#!/usr/bin/env python3
import requests
import json
import os
from typing import Dict, List
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from utils import *
import gpcrdb
from Bio.Align import PairwiseAligner
import time
from vcf import SNV
from config import AM_FILENAME, ALIGNMENT_CANDIDATES_FILENAME
import alphamissense

class BlankLineError(Exception):
    pass

class EnsemblGeneEntry:
    def __init__(self, gpcrdb_entry: gpcrdb.GPCRdbEntry, force=False):
        self.gpcrdb_entry = gpcrdb_entry

        self.ensembl_path = self.gpcrdb_entry.ensembl_path
        self.sequence_path = self.gpcrdb_entry.sequence_path
        self.visualization_path = os.path.join(self.gpcrdb_entry.dirpath, "gene.png")
        self.cds_path = self.gpcrdb_entry.cds_path
        self.alignment_path = self.gpcrdb_entry.alignment_path
        self.alphamissense_path = os.path.join(self.gpcrdb_entry.dirpath, AM_FILENAME)

        self._get_ensembl_gene_entry(force=force)
        with open(self.ensembl_path) as f:
            self._entry = json.load(f)

            if self._entry['assembly_name'] != 'GRCh38':
                raise Exception
            
            self.region = Region(normalized_chromosome(self._entry['seq_region_name']), self._entry['start'], self._entry['end'])
            self.strand = self._entry['strand']

            self.canonical_transcript_id = None
            self.canonical_translation_id = None
            for t in self._entry['Transcript']:
                if t['is_canonical'] == 1:
                    self.canonical_transcript_id = t['id']
                    self.canonical_translation_id = t['Translation']['id']
                    self.ordered_coding_regions = self._calc_coding_regions(t)
                    break
            if self.canonical_transcript_id is None:
                raise Exception("No canonical transcipt found.")
            
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

        self._assign_generic_number(force=True)
        self.generic_numbers = []
        self.segments = []
        with open(self.alignment_path) as f:
            for l in f.readlines():
                if l.startswith('#'):
                    continue
                res = MatchedResidue.from_csv_line(l)
                self.generic_numbers.append(res.generic_number)
                self.segments.append(res.segment)

    def _get_ensembl_gene_entry(self, force=False, max_try=5):
        if os.path.exists(self.ensembl_path) and force is False:
            return
        
        num_candidates = len(self.gpcrdb_entry.ensembl_id_candidates)
        if num_candidates > 1:
            print("UniProt united {} EnsEMBL IDs to this entry...".format(num_candidates))
        
        to_dump = None
        for gene_id in self.gpcrdb_entry.ensembl_id_candidates:
            uri = "https://rest.ensembl.org/lookup/id/{}?expand=1".format(gene_id)

            r = requests.get(uri, headers={"Content-Type": "application/json"})

            if not r.ok:
                try_count = 1
                retry_waits = [None, 10, 30, 60, 120, 300]
                while try_count <= max_try:
                    print("Request to EnsEMBL lookup API failed (try = {}). Waiting for {} seconds...".format(try_count, retry_waits[try_count]))
                    time.sleep(retry_waits[try_count])
                    r = requests.get(uri, headers={"Content-Type": "application/json"})
                    if r.ok:
                        break
                    try_count += 1
                else:
                    raise Exception
            
            j = r.json()

            # Check chromosome name to determine which entry to use
            try:
                _ = normalized_chromosome(j['seq_region_name'])
                if to_dump is None:
                    to_dump = j
                else:
                    raise Exception("Could not determine which entry to use...")
            except:
                print("Gene {} on {} is skipped".format(gene_id, j['seq_region_name']))
                continue
        with open(self.ensembl_path, 'w') as f:
            json.dump(to_dump, f, indent=2)

    def _get_reference_sequence(self, force=False):
        if os.path.exists(self.sequence_path) and force is False:
            return
        
        uri = "https://rest.ensembl.org/sequence/region/human"

        data = {
            "coord_system_version": "GRCh38",
            "regions": ["{}:{}..{}:{}".format(self.region.chromosome, self.region.start, self.region.end, self.strand)]
        }
        r = requests.post(uri, headers={"Content-Type": "application/json"}, data=json.dumps(data))

        if not r.ok:
            r.raise_for_status()
            raise Exception
        
        j = r.json()

        with open(self.sequence_path, 'w') as f:
            json.dump(j, f, indent=2)

    def _calc_coding_regions(self, transcript: Dict) -> List[Region]:
        """
        We only have positions where coding regions start and end
                           and where exons          start and end.
        So we must assembly exact coding regions (indicated as '*' in the following examples).
        Note that start (s) is smaller than end (e), reagrdless of strand (plus or minus).
        
                          cds_s                                cds_e
                            |                                    |
        ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
                            *****       *********       **********
         |---|     |------------|       |-------|       |--------------|        |---| 
         exon1   exon2_s      exon2_e exon3_s exon3_e exon4_s        exon4_e    exon5
        
        Or,
                          cds_s                                 cds_e
                            |                                     |
        ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
                            ***************************************
                    |------------------------------------------------------------|
                  exon6_s                                                      exon6_e
        """

        cds_s = transcript['Translation']['start']
        cds_e = transcript['Translation']['end']
        assert(cds_s < cds_e)

        coding_regions = []
        for exon in transcript['Exon']:
            exon_s = int(exon['start'])
            exon_e = int(exon['end'])
            assert(exon_s < exon_e)
            
            if exon_s <= cds_s <= exon_e <= cds_e:   # Exon 2 in the above example
                coding_regions.append(Region(self.region.chromosome, cds_s, exon_e))
            elif cds_s <= exon_s <= exon_e <= cds_e: # Exon 3
                coding_regions.append(Region(self.region.chromosome, exon_s, exon_e))
            elif cds_s <= exon_s <= cds_e <= exon_e: # Exon 4
                coding_regions.append(Region(self.region.chromosome, exon_s, cds_e))
            elif exon_s <= cds_s <= cds_e <= exon_e: # Exon 6
                coding_regions.append(Region(self.region.chromosome, cds_s, cds_e))
            else:                                    # Exons 1 and 5 contain no coding region
                continue

        coding_regions.sort(key=lambda r: r.start)

        # Remove stop codon from coding regions
        if self.strand == 1:
            assert(len(coding_regions[-1]) > 3)
            last_region = coding_regions.pop(-1)
            coding_regions.append(Region(self.region.chromosome, last_region.start, last_region.end - 3))
        else:
            assert(len(coding_regions[0]) > 3)
            last_region = coding_regions.pop(0)
            coding_regions.insert(0, Region(self.region.chromosome, last_region.start + 3, last_region.end))

        return coding_regions
    
    def _save_cds(self, force=False):
        if os.path.exists(self.cds_path) and force is False:
            return

        spliced = self._extract_spliced_coding_sequence(self.plus_strand)
        stranded = spliced if self.strand == 1 else complementary_sequence(spliced)[::-1]
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
            
            json.dump(d, f, indent=2)
            assert(stranded[:3] == 'ATG' and translated[0] == 'M')

    def _extract_spliced_coding_sequence(self, unspliced_plus_strand: str) -> str:
        offset = -self.region.start
        spliced_plus_strand = ''
        for region in self.ordered_coding_regions:
            spliced_plus_strand += unspliced_plus_strand[region.start + offset: region.end + offset + 1]
        return spliced_plus_strand

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

    def _align_residues(self):
        print("GPCRdb sequence with {} residues and Ensembl sequence with {} residues will be aligned.".format(len(self.gpcrdb_entry.protein_seq), len(self.protein_seq)))
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.1

        alignments= aligner.align(self.gpcrdb_entry.protein_seq, self.protein_seq) # target, query
        if len(alignments) > 1:
            print("Total {} alignments presented. The first one will be used.".format(len(alignments)))
            with open(os.path.join(self.gpcrdb_entry.dirpath, ALIGNMENT_CANDIDATES_FILENAME), 'w') as f:
                for alignment in alignments:
                    print(alignment, file=f)
                    print('-' * 80, file=f)
        alignment = alignments[0]

        # Assign generic numbers based on alignment
        residues = []
        for target, query in zip(*alignment.aligned):
            for gpcrdb_res_idx, ensembl_res_idx in zip(range(*target), range(*query)):
                residue = MatchedResidue(
                    self.protein_seq[ensembl_res_idx], ensembl_res_idx + 1, 
                    self.gpcrdb_entry.segments[gpcrdb_res_idx], 
                    self.gpcrdb_entry.generic_numbers[gpcrdb_res_idx], 
                    self.gpcrdb_entry.protein_seq[gpcrdb_res_idx], gpcrdb_res_idx + 1
                )
                residues.append(residue)
        print("Total {}/{} residues were filled by alignment.".format(len(residues), len(self.protein_seq)))

        # Add unaligned residues
        aligned_residue_numbers = [r.ensembl_residue_number for r in residues]
        for res_num in range(1, len(self.protein_seq) + 1):
            if res_num not in aligned_residue_numbers:
                residue = MatchedResidue(self.protein_seq[res_num - 1], res_num, None, None, None, None)
                residues.append(residue)
        residues.sort(key=lambda r: r.ensembl_residue_number)
        print("Increased from {} to {} residues by adding unaligned residues".format(len(aligned_residue_numbers), len(residues)))

        # Fill missing segments with the sandwiched segments
        for i in range(len(residues)):
            if residues[i].segment is None:
                preceding_segment = Segment.Nterm
                for residue in reversed(residues[:i]):
                    if residue.segment is not None:
                        preceding_segment = residue.segment
                        break

                following_segment = Segment.Cterm
                for residue in residues[i + 1:]:
                    if residue.segment is not None:
                        following_segment = residue.segment
                        break
                
                if preceding_segment == following_segment:
                    residues[i].segment = preceding_segment
                else:
                    # Failed to guess segment
                    residues[i].segment = Segment.FailedToGuess
        return residues

    def _assign_generic_number(self, force=False):
        if self.gpcrdb_entry.entry_name in ["ccr2_human", "gp142_human", "agrd2_human", "agrf3_human", "agrg6_human", "agrl2_human", "agrl3_human", "celr1_human"]:
            force = True

        if os.path.exists(self.alignment_path) and force is False:
            return
        
        residues = []
        if len(self.gpcrdb_entry.protein_seq) == len(self.protein_seq):
            if self.gpcrdb_entry.protein_seq != self.protein_seq:
                # In case GPCRdb and Ensembl sequences unmatched but have the same length.
                # We assume that these residues are reflecting worldwide dominant variations.
                # Strictly, this can be a combination of a insertion and a deletion with the same length.
                # However, in such cases, the number of unmatched residues will be too big.
                unmatched_residues = [1 if r1 != r2 else 0 for r1, r2 in zip(self.gpcrdb_entry.protein_seq, self.protein_seq)]
                print("CAUTION: GPCRdb and Ensembl sequence had the same length with total {} unmathced residues.".format(sum(unmatched_residues)))
            for i, (gpcrdb_aa, ensembl_aa) in enumerate(zip(self.gpcrdb_entry.protein_seq, self.protein_seq)):
                res = MatchedResidue(ensembl_aa, i + 1, 
                                     self.gpcrdb_entry.segments[i], self.gpcrdb_entry.generic_numbers[i],
                                     gpcrdb_aa, i + 1)
                residues.append(res)
        else:
            # In case insertions/deletions, alignment needed.
            residues = self._align_residues()

        with open(self.alignment_path, 'w') as f:
            lines = [MatchedResidue.header()]
            lines += [r.to_csv_line() for r in residues]
            f.write('\n'.join(lines))

    def annotate(self, snv: SNV):
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

            var_type, pathogenicity = None, None
            if alt_aa == '*':
                var_type = VariationType.NONSENSE
            elif ref_aa == alt_aa:
                var_type = VariationType.SILENT
            else:
                var_type = VariationType.MISSENSE
                substitution = ref_aa + str(res_num) + alt_aa
                try:
                    pathogenicity = alphamissense.lookup(self.alphamissense_path, snv, self.gpcrdb_entry.accession, self.canonical_transcript_id, substitution)
                except: # Could not find AlphaMissense pathogenicity due to the transcript ID mismatch
                    pathogenicity = -1
            return Annotation(var_type, snv, ref_codon, alt_codon, res_num, segment, gen_num, pathogenicity)
        raise Exception("Annotation is unavailable!")

class MatchedResidue:
    def __init__(self, ensembl_residue: str, ensembl_residue_number: int, segment: Segment, generic_number: str, gpcrdb_residue: str, gpcrdb_residue_number: int) -> None:
        self.ensembl_residue = ensembl_residue
        self.ensembl_residue_number = ensembl_residue_number
        self.segment = segment
        self.generic_number = generic_number
        self.gpcrdb_residue = gpcrdb_residue
        self.gpcrdb_residue_number = gpcrdb_residue_number

    def to_csv_line(self) -> str:
        ensembl = '-'
        gpcrdb = '-'
        if self.ensembl_residue and self.ensembl_residue_number:
            ensembl = "{}{}".format(self.ensembl_residue, self.ensembl_residue_number)
        if self.gpcrdb_residue and self.gpcrdb_residue_number:
            gpcrdb = "{}{}".format(self.gpcrdb_residue, self.gpcrdb_residue_number)
        return "{},{},{},{}".format(ensembl, self.segment, self.generic_number, gpcrdb)

    @classmethod
    def from_csv_line(cls, line: str):
        l = line.strip('\n')
        if len(l) == 0 or l.startswith('#'):
            raise BlankLineError
        cols = l.split(',')
        ensembl_residue, ensembl_residue_number = (None, None) if cols[0] == '-' else (cols[0][0], int(cols[0][1:]))
        gpcrdb_residue, gpcrdb_residue_number = (None, None) if cols[3] == '-' else (cols[3][0], int(cols[3][1:]))

        return cls(ensembl_residue, ensembl_residue_number, Segment.value_of(cols[1]), cols[2], gpcrdb_residue, gpcrdb_residue_number)
    
    @classmethod
    def header(cls):
        return "#EnsEMBL residue,Segment,Generic number,GPCRdb residue"

class Annotation:
    def __init__(self, var_type: VariationType, snv: SNV, ref_codon: str, alt_codon: str, residue_number: int, segment: Segment, generic_number: str, pathogenicity: float) -> None:
        self.var_type = var_type
        self.snv = snv

        self.ref_codon = ref_codon
        self.alt_codon = alt_codon
        assert(len(ref_codon) == len(alt_codon) == 3)
        self.ref_aa = translate(ref_codon)
        self.alt_aa = translate(alt_codon)
        self.residue_number = residue_number

        self.segment = segment
        self.generic_number = generic_number
        self.pathogenicity = pathogenicity

    def to_csv_line(self) -> str:
        cols = [self.var_type, normalized_chromosome(self.snv.chromosome), self.snv.position, self.snv.rsid]
        cols += [self.residue_number]
        cols += [self.snv.ref, self.ref_codon, self.ref_aa]
        cols += [self.snv.alt, self.alt_codon, self.alt_aa]
        cols += [self.segment, self.generic_number]
        cols += [self.snv.AC, self.snv.AN, self.snv.AF, self.pathogenicity]
        return ','.join([str(v) for v in cols])
    
    @classmethod
    def from_csv_line(cls, line: str):
        l = line.strip('\n')
        if len(l) == 0 or l.startswith('#'):
            raise BlankLineError
        cols = line.strip().split(',')
        snv = SNV(normalized_chromosome(cols[1]), int(cols[2]), cols[3], cols[5], cols[8], 
                  int(cols[13]), int(cols[14]), float(cols[15]))
        
        pathogenicity = None if cols[16] == 'None' else float(cols[16])
        generic_number = None if cols[12] == 'None' else cols[12]

        return cls(VariationType.name_of(cols[0]), snv, cols[6], cols[9], int(cols[4]),
                   Segment.value_of(cols[11]), generic_number, pathogenicity)

    @classmethod
    def header(cls) -> str:
        header = ["#Var_Type", "Chr", "Pos", "rsID", "Res_Num"]
        header += ["Ref_Base", "Ref_Codon", "Ref_AA"]
        header += ["Alt_Base", "Alt_Codon", "Alt_AA"]
        header += ["Seg", "Generic_Num"]
        header += ["AC", "AN", "AF", "pathogenicity"]
        return ','.join(header)