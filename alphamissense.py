#!/usr/bin/env python3
import gzip
import os
from typing import List
from utils import normalized_chromosome

__AM = None

class AMissense:
    def __init__(self, chromosome: str, position: int, ref: str, alt: str, variant: str, pathogenicity: float) -> None:
        assert(len(ref) == 1 and ref in "ACGT")
        assert(len(alt) == 1 and alt in "ACGT")
        assert(0 <= pathogenicity <= 1)
        self.chromosome = normalized_chromosome(chromosome)
        self.position = position
        self.ref = ref
        self.alt = alt
        self.variant = variant
        self.pathogenicity = pathogenicity

    def __str__(self) -> str:
        return "Chr{} {}, {}->{}, {}, {}".format(self.chromosome, self.position, self.ref, self.alt, self.variant, self.pathogenicity)


def extract(uniprot_id: str, fpath, force=False) -> List[AMissense]:
    if force or not os.path.exists(fpath):
        if not __AM:
            with gzip.open("./data/AlphaMissense_hg38.tsv.gz", 'rt') as f:
                __AM = f.readlines()

        with open(fpath, 'w') as f:
            for l in __AM:
                if l.startswith('#CHROM'):
                    f.write(l)
                    continue
                if uniprot_id in l:
                    f.write(l)
    
    ams = []
    with open(fpath) as f:
        for l in f:
            if l.startswith("#"):
                continue
            cols = l.strip('\n').split('\t')
            assert(cols[4] == 'hg38')
            am = AMissense(cols[0], int(cols[1]), cols[2], cols[3], cols[7], float(cols[8]))
            ams.append(am)
    return ams

if __name__ == '__main__':
    extract("O00590", "ACKR2.tsv")