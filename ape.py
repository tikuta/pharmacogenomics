import os
import json
import enum
from utils import GET_json_with_retries, Region, calc_coding_regions
import gpcrdb

class GreatApe(enum.Enum):
    # human = "homo_sapiens"
    chimpanzee = "pan_troglodytes"
    # bonobo = "pan_paniscus"
    # gorilla = "gorilla_gorilla"
    # sumatran_orangutan = "pongo_abelii"
    # bornean_orangutan = "pongo_pygmaeus" # No entries in EnsEMBL database


class EnsemblHomology:
    def __init__(self, gpcrdb_entry: gpcrdb.GPCRdbEntry, force=False) -> None:
        self.gpcrdb_entry = gpcrdb_entry
        self.dirnames = {ape: os.path.join("apes", ape.value, gpcrdb_entry.entry_name) for ape in GreatApe}
        for d in self.dirnames.values():
            os.makedirs(d, exist_ok=True)
        self._get_homologs(force=force)

        protein_ids = {ape: [] for ape in GreatApe}
        gene_ids = {ape: [] for ape in GreatApe}
        for ape in GreatApe:
            for gene_id in self.gpcrdb_entry.ensembl_id_candidates:
                p = os.path.join(self.dirnames[ape], f"{gene_id}.json")
                with open(p) as f:
                    j = json.load(f)
                    
                    if len(j['data']) == 0:
                        continue

                    assert(len(j['data']) == 1)

                    for h in j['data'][0]['homologies']:
                        protein_ids[ape].append(h['target']['protein_id'])
                        gene_ids[ape].append(h['target']['id'])
        self.protein_ids = protein_ids
        self.gene_ids = gene_ids

    def _get_homologs(self, force=False):
        for ape in GreatApe:
            for gene_id in self.gpcrdb_entry.ensembl_id_candidates:
                p = os.path.join(self.dirnames[ape], f"{gene_id}.json")
                if os.path.exists(p) and force is False:
                    continue

                uri = f'https://rest.ensembl.org/homology/id/human/{gene_id}?target_species={ape.value}'
                j = GET_json_with_retries(uri)
            
                with open(p, 'w') as f:
                    json.dump(j, f, indent=2)

class EnsemblChimpanzeeGene:
    def __init__(self, gene_id: str, force=False):
        self.gene_id = gene_id
        self.dirpath = os.path.join("apes", GreatApe.chimpanzee.value)
        self.ensembl_path = os.path.join(self.dirpath, f"{gene_id}.json")

        self._get_ensembl_gene_entry(force=force)
        with open(self.ensembl_path) as f:
            self._entry = json.load(f)

            assert(self._entry['assembly_name'] == 'Pan_tro_3.0')
            
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
        """
        self._get_reference_sequence(force=force)
        with open(self.sequence_path) as f:
            j = json.load(f)
            assert(len(j) == 1)
            if self.strand == 1:
                self.plus_strand = j[0]['seq']
            else:
                self.plus_strand = complementary_sequence(j[0]['seq'])[::-1]

        self._save_cds(force=force)
        with open(self.cds_path) as f:
            j = json.load(f)
            self.protein_seq = j['translated']
            self.stranded_coding_sequence = j['canonical_cds']

        self._assign_generic_number(force=force)
        self.generic_numbers = []
        self.segments = []
        with open(self.alignment_path) as f:
            for l in f.readlines():
                if l.startswith('#'):
                    continue
                res = MatchedResidue.from_csv_line(l)
                self.generic_numbers.append(res.generic_number)
                self.segments.append(res.segment)
        """
    def _get_ensembl_gene_entry(self, force=False):
        if os.path.exists(self.ensembl_path) and force is False:
            return
        
        uri = f"https://rest.ensembl.org/lookup/id/{self.gene_id}?expand=1"
        j = GET_json_with_retries(uri)

        with open(self.ensembl_path, 'w') as f:
            json.dump(j, f, indent=2)

    def _assign_generic_number(self, force=False):
        pass