import os
import json
import enum
from utils import *
import gpcrdb
from config import REFERENCE_CHIMPANZEE_GENOME
from dataclasses import dataclass
from typing import List
import vcf
import genbank

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

@dataclass
class Annotation:
    var_type: VariationType
    snv: vcf.SNV

    ref_codon: str
    alt_codon: str

    residue_number: int
    segment: Segment
    generic_number: str

    def __post_init__(self):
        assert(len(self.ref_codon) == len(self.alt_codon) == 3)

        self.ref_aa = translate(self.ref_codon)
        self.alt_aa = translate(self.alt_codon)

    def to_csv_line(self) -> str:
        cols = [self.var_type, self.snv.chromosome, self.snv.position, self.snv.rsid]
        cols += [self.residue_number]
        cols += [self.snv.ref, self.ref_codon, self.ref_aa]
        cols += [self.snv.alt, self.alt_codon, self.alt_aa]
        cols += [self.segment, self.generic_number]
        cols += [self.snv.AC, self.snv.AN, self.snv.AF]
        return ','.join([str(v) for v in cols])
    
    @classmethod
    def from_csv_line(cls, line: str):
        l = line.strip('\n')
        if len(l) == 0 or l.startswith('#'):
            return None
        cols = line.strip().split(',')
        snv = vcf.SNV(cols[1], int(cols[2]), cols[3], cols[5], cols[8], 
                      int(cols[13]), int(cols[14]), float(cols[15]))
        
        generic_number = None if cols[12] == 'None' else cols[12]

        return cls(VariationType.name_of(cols[0]), snv, cols[6], cols[9], int(cols[4]),
                   Segment.value_of(cols[11]), generic_number)

    @classmethod
    def header(cls) -> str:
        header = ["#Var_Type", "Chr", "Pos", "rsID", "Res_Num"]
        header += ["Ref_Base", "Ref_Codon", "Ref_AA"]
        header += ["Alt_Base", "Alt_Codon", "Alt_AA"]
        header += ["Seg", "Generic_Num"]
        header += ["AC", "AN", "AF"]
        return ','.join(header)

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
        self.sequence_path = os.path.join(self.dirpath, f"{self.gene_id}_seq.json")
        self.annotated_csv_path = os.path.join(self.dirpath, f"{self.gene_id}_CDS.csv")
        self.cds_path = os.path.join(self.dirpath, f"{self.gene_id}_CDS.json")
        self.genbank_path = os.path.join(self.dirpath, f"{self.gene_id}.gb")
        self.visualization_path = os.path.join(self.dirpath, f"{self.gene_id}.png")

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
                    if self.gene_id == 'ENSPTRG00000002704': # OPN4 has illegal reading frame
                        self.ordered_coding_regions[-1].end += 2
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

        self._get_reference_sequence(force=force)
        with open(self.sequence_path) as f:
            j = json.load(f)
            assert(len(j) == 1)
            if self.strand == 1:
                self.plus_strand = j[0]['seq']
            else:
                self.plus_strand = complementary_sequence(j[0]['seq'])[::-1]

        self._save_genbank(force=force)
        self.visualize()
        self._save_cds(force=force)
        with open(self.cds_path) as f:
            j = json.load(f)
            self.protein_seq = j['translated']
            self.stranded_coding_sequence = j['canonical_cds']

        self._assign_generic_number()
        self.generic_numbers = []
        self.segments = []
        with open(self.alignment_path) as f:
            for l in f.readlines():
                res = MatchedResidue.from_csv_line(l)
                if res:
                    self.generic_numbers.append(res.generic_number)
                    self.segments.append(res.segment)

    def visualize(self, force=False):
        if os.path.exists(self.visualization_path) and force is False:
            return

        fig, ax = plt.subplots(1, 1, figsize=(2, 5), dpi=150)

        ax.set_yticks([self.region.start, self.region.end])
        ax.set_yticklabels([self.region.start, self.region.end], size=6)

        ticks = []
        total_bp = 0
        for region in self.ordered_coding_regions:
            ax.axhspan(region.start, region.end, color='tab:orange', zorder=100)
            bp = region.end - region.start + 1
            total_bp += bp
            text = "{} bp".format(bp)
            ax.text(0, (region.start + region.end) / 2, text, size=6, va='center', ha='center', zorder=1000)
            ticks.append(region.start)
            ticks.append(region.end)

        ax.set_xlim(-0.4, 0.4)
        ax.set_ylim(self.region.start, self.region.end)
        ax.set_yticks(ticks, minor=True)
        ax.set_yticklabels(ticks, size=6, minor=True)

        label = "Total {} bp / {} aa".format(total_bp, total_bp / 3)
        ax.bar(0, self.region.end, 0.8, tick_label=label, color='tab:gray', zorder=10)

        ax.set_title("chr{} ({})".format(self.region.chromosome, "+" if self.strand == 1 else "-"), size=8)

        fig.tight_layout()

        fig.savefig(self.visualization_path)
        plt.close(fig)

    def _get_ensembl_gene_entry(self, force=False):
        if os.path.exists(self.ensembl_path) and force is False:
            return
        
        uri = f"https://rest.ensembl.org/lookup/id/{self.gene_id}?expand=1"
        j = GET_json_with_retries(uri)

        with open(self.ensembl_path, 'w') as f:
            json.dump(j, f, indent=2)

    def _get_reference_sequence(self, force=False):
        if self.species != GreatApe.chimpanzee:
            raise NotImplementedError

        if os.path.exists(self.sequence_path) and force is False:
            return
        
        uri = f"https://rest.ensembl.org/sequence/region/{self.species.value}"

        data = {
            "coord_system_version": REFERENCE_CHIMPANZEE_GENOME,
            "regions": ["{}:{}..{}:{}".format(self.region.chromosome, self.region.start, self.region.end, self.strand)]
        }

        j = POST_with_data_to_get_json_with_retries(uri, data)

        with open(self.sequence_path, 'w') as f:
            json.dump(j, f, indent=2)

    def _save_genbank(self, force=False):
        force = True
        if os.path.exists(self.genbank_path) and force is False:
            return
        
        labels = []
        is_forward = True if self.strand == 1 else False
        for r in self.ordered_coding_regions:
            name = f"CDS ({self.region.chromosome}:{r.start}-{r.end})"
            label = genbank.GenBankLabel(name, r.start - self.region.start + 1, r.end - self.region.start + 1, is_forward)
            labels.append(label)
        gb = genbank.GenBank(self.plus_strand, labels)
        gb.save(self.genbank_path, overwrite=True)

    def _save_cds(self, force=False):
        if os.path.exists(self.cds_path) and force is False:
            return

        spliced = self._extract_spliced_coding_sequence(self.plus_strand)
        stranded = spliced if self.strand == 1 else complementary_sequence(spliced)[::-1]
        if len(stranded) % 3 != 0:
            for r in [self.region] + self.ordered_coding_regions:
                print(r.start, r.end, r.end - r.start + 1, "bp")
            gb = genbank.GenBank(stranded, [])
            gb.save("error.gb", overwrite=True)
            raise Exception(f"{self.gene_id}", len(stranded), spliced)
        translated = translate(stranded)
        with open(self.cds_path, 'w') as f:
            d = {
                "gene_region": "{}:{}".format(str(self.region), self.strand),
                "canonical_cds_region": [str(r) for r in self.ordered_coding_regions],
                "canonical_transcript_id": self.canonical_transcript_id,
                "canonical_translation_id": self.canonical_translation_id,
                "plus_strand": spliced,
                "canonical_cds": stranded,
                "translated": translated
            }
            
            if stranded[:3] != 'ATG' or translated[0] != 'M':
                raise Exception(f"{self.gene_id}")
            json.dump(d, f, indent=2)

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
    
    def _extract_spliced_coding_sequence(self, unspliced_plus_strand: str) -> str:
        offset = -self.region.start
        spliced_plus_strand = ''
        for region in self.ordered_coding_regions:
            spliced_plus_strand += unspliced_plus_strand[region.start + offset: region.end + offset + 1]
        return spliced_plus_strand

    def annotate(self, snv: vcf.SNV) -> Annotation:
        assert(self.region.chromosome == snv.chromosome)
        assert(self.region.start <= snv.position <= self.region.end)
        assert(True in [r.start <= snv.position <= r.end for r in self.ordered_coding_regions])

        offset = -self.region.start
        altered_plus_strand = self.plus_strand[:snv.position + offset] + snv.alt + self.plus_strand[snv.position + offset + 1:]
        altered_spliced = self._extract_spliced_coding_sequence(altered_plus_strand)
        altered_stranded = altered_spliced if self.strand == 1 else complementary_sequence(altered_spliced)[::-1]
        altered_translated = translate(altered_stranded)
        
        for res_num in range(1, len(self.protein_seq) + 1):
            ref_codon = self.stranded_coding_sequence[3 * (res_num - 1): 3 * res_num]
            alt_codon = altered_stranded[3 * (res_num - 1): 3 * res_num]
            ref_aa = self.protein_seq[res_num - 1]
            alt_aa = altered_translated[res_num - 1]
            gen_num = self.generic_numbers[res_num - 1]
            segment = self.segments[res_num - 1]

            if ref_codon == alt_codon:
                continue

            var_type = None
            if alt_aa == '*':
                var_type = VariationType.NONSENSE
            elif ref_aa == alt_aa:
                var_type = VariationType.SILENT
            else:
                var_type = VariationType.MISSENSE
            return Annotation(var_type, snv, ref_codon, alt_codon, res_num, segment, gen_num)
        raise Exception("Annotation is unavailable!")
