#!/usr/bin/env python3
from typing import List
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('Agg')
import os
import enum
from misc import STANDARD_CODES
import json

class Segment(enum.Enum):
    Nterm = "N-term"
    TM1 = "TM1"
    ICL1 = "ICL1"
    TM2 = "TM2"
    ECL1 = "ECL1"
    TM3 = "TM3"
    ICL2 = "ICL2"
    TM4 = "TM4"
    ECL2 = "ECL2"
    TM5 = "TM5"
    ICL3 = "ICL3"
    TM6 = "TM6"
    ECL3 = "ECL3"
    TM7 = "TM7"
    ICL4 = "ICL4"
    H8 = "H8"
    Cterm = "C-term"
    NONE = "None"

    def __str__(self) -> str:
        return self.value
    
    @property
    def color(self):
        cmap = plt.get_cmap('rainbow', len(Segment) - 1) 
        if self == Segment.NONE:
            return 'tab:gray'
        for i, seg in enumerate(Segment):
            if seg == self:
                return cmap(i)
            
    @property
    def index(self):
        for i, seg in enumerate(Segment):
            if seg == self:
                return i
    
    @classmethod
    def value_of(cls, target):
        for e in cls:
            if e.value == target:
                return e
        raise ValueError

    @classmethod
    def generic_number_of(cls, g_num):
        seg = int(g_num.split('.')[0])
        if 1 <= seg <= 7:
            return cls.value_of("TM" + str(seg))
        elif seg == 8:
            return cls.H8
        elif seg == 12:
            return cls.ICL1
        elif seg == 23:
            return cls.ECL1
        elif seg == 34:
            return cls.ICL2
        elif seg == 45:
            return cls.ECL2
        elif seg == 56:
            return cls.ICL3
        elif seg == 67:
            return cls.ECL3
        raise ValueError

class VariationType(enum.Enum):
    MISSENSE = 1
    SILENT = 0
    NONSENSE = -1

    def __str__(self) -> str:
        return self.name
    
    @classmethod
    def name_of(cls, target):
        for e in cls:
            if e.name == target:
                return e
        raise ValueError
    
class Annotation:
    def __init__(self, ref_codon: str, alt_codon: str, res_num: int) -> None:
        assert(len(ref_codon) == len(alt_codon) == 3)
        self.ref_codon = ref_codon
        self.ref_aa = translate(ref_codon)
        self.alt_codon = alt_codon
        self.alt_aa = translate(alt_codon)
        self.res_num = res_num
        if self.ref_aa == self.alt_aa:
            self.var_type = VariationType.SILENT
        else:
            if self.alt_aa == '*':
                self.var_type = VariationType.NONSENSE
            else:
                self.var_type = VariationType.MISSENSE

class Region:
    def __init__(self, chromosome: str, start: int, end: int) -> None:
        self.chromosome = normalized_chromosome(chromosome)
        assert(start <= end)
        self.start = start
        self.end = end

    def __str__(self) -> str:
        return "{}:{}-{}".format(self.chromosome, self.start, self.end)
    
    def __len__(self) -> int:
        return self.end - self.start + 1

def normalized_chromosome(c) -> str:
    if isinstance(c, int):
        return str(c)
    
    if c.startswith('chr'):
        return c[3:]
    
    if c not in [str(v) for v in range(1, 23)] + ['X', 'Y']:
        print("Chromosome `{}` is not in 38KJPN calls".format(c))
    
    return c

def complementary_sequence(seq: str) -> str:
    pairs = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return ''.join([pairs[s] for s in seq])

def translate(seq: str, **unusual_codons) -> str:
    assert(len(seq) % 3 == 0)
    
    CODES = {**STANDARD_CODES, **unusual_codons}

    ret = ''
    for i in range(0, len(seq), 3):
        triplet = seq[i: i + 3]
        ret += CODES[triplet]
    return ret

class SNV:
    def __init__(self, chromosome: str, position: int, rsid: str, ref: str, alt: str, AC, AN, AF) -> None:
        assert(len(ref) == len(alt) == 1)
        self.chromosome = normalized_chromosome(chromosome)
        self.position = position
        self.rsid = rsid
        self.ref = ref
        self.alt = alt

        self.AC = AC
        self.AN = AN
        self.AF = AF

class Variation:
    def __init__(self, line: str) -> None:
        cols = line.strip().split('\t')

        self.chromosome = normalized_chromosome(cols[0])
        self.position = int(cols[1])
        self.rsid = cols[2]
        self.ref = cols[3]
        assert(',' not in self.ref)
        self.alts = cols[4].split(',')
        try:
            self.quality = float(cols[5])
        except:
            self.quality = cols[5]
        self.passed = True if cols[6] == 'PASS' else False
        
        # Some keys may not have values (e.g., TOMMO_POSSIBLE_PLATFORM_BIAS_SITE)
        keyvals = {keyval.split('=')[0]: keyval.split('=')[-1] for keyval in cols[7].split(';')}

        self.AC = [int(v) for v in keyvals['AC'].split(',')]
        self.AN = int(keyvals['AN'])
        self.AF = keyvals['AF'].split(',') # retain as string

        assert(len(self.alts) == len(self.AC) == len(self.AF))

    def snvs(self) -> List[SNV]:
        ret = []
        for i, alt in enumerate(self.alts):
            if len(self.ref) == 1:
                assert(self.ref in 'ACGT')
                if len(alt) == 1: # 1:1 SNV
                    assert(alt in 'ACGT') # Illegal VCF can contain '-' to indicate deletions
                    ret.append(SNV(self.chromosome, self.position, self.rsid, self.ref, alt, self.AC[i], self.AN, self.AF[i]))
                else: # 1:N insertion
                    pass
            else:
                if len(alt) == 1: # N:1 deletion
                    pass
                elif len(self.ref) == len(alt): # N:N substitution
                    count = 0
                    idx = None
                    for j, (b1, b2) in enumerate(zip(self.ref, alt)):
                        if b1 != b2:
                            count += 1
                            idx = j
                    if count == 1: # 1:1 SNV
                        ret.append(SNV(self.chromosome, self.position + idx, self.rsid, self.ref[idx], alt[idx], self.AC[i], self.AN, self.AF[i]))
                    else:
                        # Combination of 1:1 SNV (e.g. ACGT->ACCC)
                        # No such calls in 54KJPN
                        raise NotImplementedError
                else: # N:M indel
                    pass
        return ret

class Gene:
    def __init__(self, region: Region, seq: str, strand: int, coding_regions: List[Region]) -> None:
        assert(len(seq) == region.end - region.start + 1)
        self.region = region

        coding_regions.sort(key=lambda r: r.start)
        assert(region.start <= coding_regions[0].start <= coding_regions[-1].end <= region.end)
        self.coding_regions = coding_regions

        assert(strand == 1 or strand == -1)
        self.strand = strand
        if strand == 1:
            self.plus = seq
        else:
            self.plus = complementary_sequence(seq)[::-1]

    def __len__(self) -> int:
        return len(self.region)
        
    def at(self, region: Region, strand: int, snv:SNV=None) -> str:
        assert(self.region.chromosome == region.chromosome)
        assert(self.region.start <= region.start <= region.end <= self.region.end)
        
        plus = self.plus
        if snv:
            assert(self.region.start <= snv.position <= self.region.end)
            idx = snv.position - self.region.start
            assert(self.at(Region(snv.chromosome, snv.position, snv.position), 1) == snv.ref)
            plus = plus[:idx] + snv.alt + plus[idx + 1:]

        offset = region.start - self.region.start
        length = region.end - region.start + 1
        if strand == 1:
            return plus[offset: offset + length]
        else:
            return complementary_sequence(plus[offset: offset + length])[::-1]
    
    def save_cds(self, fpath, force=False, snv:SNV=None):
        if os.path.exists(fpath) and force is False:
            return
        
        nucleotides = ''
        for region in self.coding_regions[::self.strand]:
            nucleotides += self.at(region, self.strand, snv=snv)
        with open(fpath, 'w') as f:
            json.dump({"nucleotides": nucleotides}, f)

    def translate(self, snv:SNV=None) -> str:
        t = ''

        for region in self.coding_regions[::self.strand]:
            t += self.at(region, self.strand, snv=snv)

        return translate(t)
    
    def annotate(self, snv: SNV) -> Annotation:
        if self.region.chromosome != snv.chromosome:
            print(self.region.chromosome, snv.chromosome)
        assert(self.region.chromosome == snv.chromosome)
        assert(self.region.start <= snv.position <= self.region.end)

        ref_t = ''
        alt_t = ''
        for region in self.coding_regions[::self.strand]:
            ref = self.at(region, self.strand)
            ref_t += ref
            if region.start <= snv.position <= region.end:
                alt_t += self.at(region, self.strand, snv=snv)
            else:
                alt_t += ref

        assert(len(ref_t) == len(alt_t) and len(ref_t) % 3 == 0)

        total_res_num = int(len(ref_t) / 3)
        for res_num in range(1, total_res_num + 1):
            ref_codon = ref_t[3 * (res_num - 1): 3 * res_num]
            alt_codon = alt_t[3 * (res_num - 1): 3 * res_num]

            if ref_codon != alt_codon:
                return Annotation(ref_codon, alt_codon, res_num)

        raise Exception
    
    def visualize(self, fpath, force=False):
        if os.path.exists(fpath) and force is False:
            return

        fig, ax = plt.subplots(1, 1, figsize=(2, 5), dpi=150)

        ax.set_yticks([self.region.start, self.region.end])
        ax.set_yticklabels([self.region.start, self.region.end], size=6)

        ticks = []
        total_bp = 0
        for region in self.coding_regions:
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

        fig.savefig(fpath)
        plt.close(fig)

if __name__ == '__main__':
    # r = Region('Z', 1, 12)
    # seq = Sequence(r, 'ACCG' * 3, -1)
    # print(seq.at(Region("Z", 1, 3), -1))
    var_type = VariationType.MISSENSE
    print(var_type.name)