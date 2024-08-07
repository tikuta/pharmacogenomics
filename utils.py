#!/usr/bin/env python3
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('Agg')
import enum
from config import STANDARD_CODES
import requests
import json
import time
from typing import List

class GproteinCoupling(enum.Enum):
    Gs = "Gs"
    Gio = "Gi/o"
    Gq11 = "Gq/11"
    G1213 = "G12/13"
    Unknown = "Unknown"

    def __str__(self) -> str:
        return self.value
    
    @property
    def color(self):
        if self == GproteinCoupling.Gs:
            return 'cyan'
        elif self == GproteinCoupling.Gio:
            return 'magenta'
        elif self == GproteinCoupling.Gq11:
            return 'yellow'
        elif self == GproteinCoupling.G1213:
            return 'white'
        elif self == GproteinCoupling.Unknown:
            return 'tab:gray'

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
    FailedToGuess = "Failed to guess"

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
    def terms(cls):
        return [cls.Nterm, cls.Cterm]
    
    @classmethod
    def TMs(cls):
        return [cls.TM1, cls.TM2, cls.TM3, cls.TM4, cls.TM5, cls.TM6, cls.TM7]
    
    @classmethod
    def helices(cls):
        return cls.TMs() + [cls.H8]

    @classmethod
    def ICLs(cls):
        return [cls.ICL1, cls.ICL2, cls.ICL3, cls.ICL4]
    
    @classmethod
    def ECLs(cls):
        return [cls.ECL1, cls.ECL2, cls.ECL3]
    
    @classmethod
    def value_of(cls, target):
        for e in cls:
            if e.value == target:
                return e
        raise ValueError

    @classmethod
    def generic_number_of(cls, g_num):
        if '.' in g_num:
            seg = int(g_num.split('.')[0])
        elif 'x' in g_num:
            seg = int(g_num.split('x')[0])
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
    def name_of(cls, target: str):
        for e in cls:
            if e.name == target:
                return e
        raise ValueError

class Region:
    def __init__(self, chromosome: str, start: int, end: int) -> None:
        self.chromosome = chromosome
        assert(start <= end)
        self.start = start
        self.end = end

    @classmethod
    def load(cls, s: str):
        elems = s.split(':')
        assert(2 <= len(elems) <= 3)

        chromosome = normalized_chromosome(elems[0])

        arange = elems[1].split('-')
        assert(len(arange) == 2)
        start, end = int(arange[0]), int(arange[1])

        strand = int(elems[2]) if len(elems) == 3 else None

        return cls(chromosome, start, end), strand

    def __str__(self) -> str:
        return "{}:{}-{}".format(self.chromosome, self.start, self.end)
    
    def __len__(self) -> int:
        return self.end - self.start + 1

def normalized_chromosome(c, species="human"):
    if species == "human":
        return normalized_human_chromosome(c)
    elif species == "chimpanzee":
        return normalized_chimpanzee_chromosome(c)
    raise Exception(f"Unknown species: `{c}`.")

def normalized_human_chromosome(c) -> str:
    if isinstance(c, int):
        return str(c)
    
    if isinstance(c, str) and c.startswith('chr'):
        return c[3:]
    
    if c in [str(v) for v in range(1, 23)] + ['X', 'Y']:
        return c

    raise Exception(f"Unknown human chromosome: `{c}`.")

def normalized_chimpanzee_chromosome(c) -> str:
    if isinstance(c, int):
        return str(c)
    
    if isinstance(c, str) and c.startswith('chr'):
        return c[3:]
    
    if [str(v) for v in range(1, 23)] + ['2A', '2B', 'X', 'Y']:
        return c

    raise Exception(f"Unknown chimpanzee chromosome: `{c}`.")

def complementary_sequence(seq: str) -> str:
    pairs = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    return ''.join([pairs[s] for s in seq])

def translate(seq: str, **unusual_codons) -> str:
    assert(len(seq) % 3 == 0)
    
    CODES = {**STANDARD_CODES, **unusual_codons}

    ret = ''
    for i in range(0, len(seq), 3):
        triplet = seq[i: i + 3]
        ret += CODES[triplet]
    return ret

def calc_coding_regions(gene_region: Region, strand: int, transcript: dict) -> List[Region]:
    """
    We only have positions where coding regions start and end
                       and where exons          start and end.
    So we must assembly exact coding regions (indicated as '*' in the following examples).
    Note that start (s) is always smaller than end (e), reagrdless of strand (plus or minus).
    
                        cds_s                                cds_e
                        |                                    |
    NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
                        *****       *********       **********
     |---|     |------------|       |-------|       |--------------|        |---| 
     exon1   exon2_s      exon2_e exon3_s exon3_e exon4_s        exon4_e    exon5
    
    Or,
                        cds_s                                 cds_e
                        |                                     |
    NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
                        ***************************************
                |------------------------------------------------------------|
                exon6_s                                                      exon6_e
    """

    cds_s = transcript['Translation']['start']
    cds_e = transcript['Translation']['end']
    assert(cds_s <= cds_e)

    coding_regions = []
    for exon in transcript['Exon']:
        exon_s = int(exon['start'])
        exon_e = int(exon['end'])
        assert(exon_s <= exon_e)
        
        if exon_s <= cds_s <= exon_e <= cds_e:   # Exon 2 in the above example
            coding_regions.append(Region(gene_region.chromosome, cds_s, exon_e))
        elif cds_s <= exon_s <= exon_e <= cds_e: # Exon 3
            coding_regions.append(Region(gene_region.chromosome, exon_s, exon_e))
        elif cds_s <= exon_s <= cds_e <= exon_e: # Exon 4
            coding_regions.append(Region(gene_region.chromosome, exon_s, cds_e))
        elif exon_s <= cds_s <= cds_e <= exon_e: # Exon 6
            coding_regions.append(Region(gene_region.chromosome, cds_s, cds_e))
        else:                                    # Exons 1 and 5 contain no coding region
            continue

    coding_regions.sort(key=lambda r: r.start)

    # Remove stop codon from coding regions
    if strand == 1:
        assert(len(coding_regions[-1]) > 3)
        last_region = coding_regions.pop(-1)
        coding_regions.append(Region(gene_region.chromosome, last_region.start, last_region.end - 3))
    else:
        assert(len(coding_regions[0]) > 3)
        last_region = coding_regions.pop(0)
        coding_regions.insert(0, Region(gene_region.chromosome, last_region.start + 3, last_region.end))

    return coding_regions

def GET_json_with_retries(uri, retry_waits=[10, 30, 60, 120, 300]):
    for count, wait in enumerate(retry_waits):
        if count > 0:
            print(f"Request to EnsEMBL lookup API ({uri}) failed (try = {count})."
                  f"Waiting for {wait} seconds...")
            time.sleep(wait)
        r = requests.get(uri, headers={"Content-Type": "application/json"})

        if r.ok:
            return r.json()
    r.raise_for_status()
        
def POST_with_data_to_get_json_with_retries(uri, data, retry_waits=[10, 30, 60, 120, 300]):
    for count, wait in enumerate(retry_waits):
        if count > 0:
            print(f"Request to EnsEMBL lookup API ({uri}) failed (try = {count})."
                  f"Waiting for {wait} seconds...")
            time.sleep(wait)
        r = requests.post(uri, headers={"Content-Type": "application/json"}, data=json.dumps(data))

        if r.ok:
            return r.json()
    r.raise_for_status()