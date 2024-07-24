import os
import json
import enum
from utils import GET_json_with_retries, Region, calc_coding_regions, Segment
import gpcrdb
from config import REFERENCE_CHIMPANZEE_GENOME
from dataclasses import dataclass
from typing import List

class GreatApe(enum.Enum):
    # human = "homo_sapiens"
    chimpanzee = "pan_troglodytes"
    # bonobo = "pan_paniscus"
    # gorilla = "gorilla_gorilla"
    # sumatran_orangutan = "pongo_abelii"
    # bornean_orangutan = "pongo_pygmaeus" # No entries in EnsEMBL database

@dataclass
class MatchedResidue:
    source_residue_aa: str
    source_residue_number: int

    target_residue_aa: str
    target_residue_number: int

    generic_number: str
    segment: Segment

    def to_csv_line(self) -> str:
        source_residue, target_residue = '-', '-'
        if self.source_residue_aa and self.source_residue_number:
            source_residue = f"{self.source_residue_aa}{self.source_residue_number}"
        if self.target_residue_aa and self.target_residue_number:
            target_residue = f"{self.target_residue_aa}{self.target_residue_number}"
        return ",".join([source_residue, str(self.segment), str(self.generic_number), target_residue])

    @classmethod
    def from_csv_line(cls, line: str):
        l = line.strip('\n')
        if len(l) == 0 or l.startswith('#'):
            return None
        cols = l.split(',')
        source_residue_aa, source_residue_number = (None, None) if cols[0] == '-' else (cols[0][0], int(cols[0][1:]))
        target_residue_aa, target_residue_number = (None, None) if cols[3] == '-' else (cols[3][0], int(cols[3][1:]))

        return cls(source_residue_aa, source_residue_number, 
                   target_residue_aa, target_residue_number, 
                   cols[2], Segment.value_of(cols[1]))
    
    @classmethod
    def header(cls):
        return "#Source residue,Segment,Generic number,Target residue"

class EnsemblHomology:
    def __init__(self, gpcrdb_entry: gpcrdb.GPCRdbEntry, species: enum.Enum, force=False) -> None:
        self.gpcrdb_entry = gpcrdb_entry
        self.species = species
        self.dirname = os.path.join("apes", species.value, gpcrdb_entry.entry_name)
        os.makedirs(self.dirname, exist_ok=True)
        self._get_homologs(force=force)

        self.gene_ids = {}
        for gene_id in self.gpcrdb_entry.ensembl_id_candidates:
            p = os.path.join(self.dirname, f"{gene_id}.json")
            with open(p) as f:
                j = json.load(f)
                
                if len(j['data']) == 0:
                    continue

                assert(len(j['data']) == 1)

                for h in j['data'][0]['homologies']:
                    self.gene_ids[gene_id] = h['target']['id']


    def _get_homologs(self, force=False):
        for gene_id in self.gpcrdb_entry.ensembl_id_candidates:
            p = os.path.join(self.dirname, f"{gene_id}.json")
            if os.path.exists(p) and force is False:
                continue

            uri = f'https://rest.ensembl.org/homology/id/human/{gene_id}?target_species={self.species.value}'
            j = GET_json_with_retries(uri)
        
            with open(p, 'w') as f:
                json.dump(j, f, indent=2)

class EnsemblHomologGene:
    def __init__(self, gene_id: str, species: enum.Enum, human_gpcrdb_entry: gpcrdb.GPCRdbEntry, human_gene_id: str, force=False):
        self.gene_id = gene_id
        self.human_gpcrdb_entry = human_gpcrdb_entry
        self.human_gene_id = human_gene_id
        self.species = species
        self.dirpath = os.path.join("apes", species.value)
        self.ensembl_path = os.path.join(self.dirpath, f"{gene_id}.json")
        self.homolog_path = os.path.join(self.dirpath, self.human_gpcrdb_entry.entry_name, f"{self.human_gene_id}.json")
        self.alignment_path = os.path.join(self.dirpath, f"{self.gene_id}.csv")

        self._get_ensembl_gene_entry(force=force)
        with open(self.ensembl_path) as f:
            self._entry = json.load(f)

            if self.species == GreatApe.chimpanzee:
                assert(self._entry['assembly_name'] == REFERENCE_CHIMPANZEE_GENOME)
            
            self.region = Region(self._entry['seq_region_name'], self._entry['start'], self._entry['end'])
            self.strand = self._entry['strand']

            self.canonical_transcript_id = None
            self.canonical_translation_id = None
            for t in self._entry['Transcript']:
                if t['is_canonical'] == 1:
                    self.canonical_transcript_id = t['id']
                    self.canonical_translation_id = t['Translation']['id']
                    self.ordered_coding_regions = calc_coding_regions(self.region, self.strand, t)
                    break
            if self.canonical_transcript_id is None:
                raise Exception("No canonical transcipt found.")
            
        with open(self.homolog_path) as f:
            j = json.load(f)
            assert(len(j['data']) == 1)

            self.homolog_align_seq, self.homolog_cigar_line = None, None
            self.human_align_seq, self.human_cigar_line = None, None

            homologies = j['data'][0]['homologies']
            for homology in homologies:
                if homology['target']['id'] == self.gene_id:
                    self.homolog_align_seq = homology['target']['align_seq']
                    self.human_align_seq = homology['source']['align_seq']

        self._assign_generic_number()

    def _get_ensembl_gene_entry(self, force=False):
        if os.path.exists(self.ensembl_path) and force is False:
            return
        
        uri = f"https://rest.ensembl.org/lookup/id/{self.gene_id}?expand=1"
        j = GET_json_with_retries(uri)

        with open(self.ensembl_path, 'w') as f:
            json.dump(j, f, indent=2)

    def _assign_generic_number(self, force=False):
        if os.path.exists(self.alignment_path) and force is False:
            return
        
        human_residues: List[MatchedResidue] = []
        with open(self.human_gpcrdb_entry.alignment_path) as f:
            for l in f:
                residue = MatchedResidue.from_csv_line(l)
                if residue and residue.source_residue_aa:
                    human_residues.append(residue)

        new_matched_residues: List[MatchedResidue] = []
        for idx, homolog_aa in enumerate(self.homolog_align_seq):
            if homolog_aa == '-':
                continue
            homolog_resnum = idx - self.homolog_align_seq[:idx].count('-') + 1

            human_aa = self.human_align_seq[idx]
            human_resnum = idx - self.human_align_seq[:idx].count('-') + 1

            if human_aa != '-':
                # look up matched residue
                for mr in human_residues:
                    if mr.source_residue_aa == human_aa and mr.source_residue_number == human_resnum:
                        new_mr = MatchedResidue(homolog_aa, homolog_resnum,
                                                human_aa, human_resnum,
                                                mr.generic_number, mr.segment)
                        new_matched_residues.append(new_mr)
                        break
            else:
                new_mr = MatchedResidue(homolog_aa, homolog_resnum, None, None, None, None)
                new_matched_residues.append(new_mr)
        
        with open(self.alignment_path, 'w') as f:
            f.write(MatchedResidue.header() + '\n')
            f.write('\n'.join([mr.to_csv_line() for mr in new_matched_residues]))