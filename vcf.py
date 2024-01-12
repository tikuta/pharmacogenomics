#!/usr/bin/env python3

import subprocess
import os
from utils import Region
from typing import List, Iterator
from utils import normalized_chromosome

class NotPassedError(Exception):
    pass

class BlankLineError(Exception):
    pass

class VariationEntry:
    def __init__(self, chromosome: str, position: int, rsid: str, ref: str, alts: List[str], ACs: List[int], AN: int, AFs: List[float]):
        self.chromosome = normalized_chromosome(chromosome)
        self.position = position
        self.rsid = rsid
        self.ref = ref
        self.alts = alts
        self.ACs = ACs
        self.AN = AN
        self.AFs = AFs
    
    @classmethod
    def load_from_54KJPN(cls, line: str):
        l = line.strip()
        if len(l) == 0:
            raise BlankLineError
        cols = l.split('\t')

        passed = True if cols[6] == 'PASS' else False
        if not passed:
            raise NotPassedError(cols[6])

        chromosome = cols[0]
        position = int(cols[1])
        rsid = cols[2]
        ref = cols[3]
        assert(',' not in ref)
        alts = cols[4].split(',')
        
        # Some keys may not have values (e.g., TOMMO_POSSIBLE_PLATFORM_BIAS_SITE)
        keyvals = {keyval.split('=')[0]: keyval.split('=')[-1] for keyval in cols[7].split(';')}

        ACs = [int(v) for v in keyvals['AC'].split(',')]
        AN = int(keyvals['AN'])
        AFs = [float(v) for v in keyvals['AF'].split(',')]

        assert(len(alts) == len(ACs) == len(AFs))

        return cls(chromosome, position, rsid, ref, alts, ACs, AN, AFs)
    
    @classmethod
    def load_from_1KGP(cls, line: str):
        l = line.strip()
        if len(l) == 0:
            raise BlankLineError
        cols = l.split('\t')[:8]

        passed = True if cols[6] == 'PASS' else False
        if not passed:
            raise NotPassedError(cols[6])

        chromosome = cols[0]
        position = int(cols[1])
        rsid = None
        ref = cols[3]
        assert(',' not in ref)
        alts = cols[4].split(',')
        
        keyvals = {keyval.split('=')[0]: keyval.split('=')[-1] for keyval in cols[7].split(';')}

        ACs = [int(v) for v in keyvals['AC'].split(',')]
        AN = int(keyvals['AN'])
        AFs = [float(v) for v in keyvals['AF'].split(',')]

        assert(len(alts) == len(ACs) == len(AFs))

        return cls(chromosome, position, rsid, ref, alts, ACs, AN, AFs)

    @property
    def snvs(self):
        for i, alt in enumerate(self.alts):
            if len(self.ref) == 1:
                assert(self.ref in 'ACGT')
                if len(alt) == 1: # 1:1 SNV
                    assert(alt in 'ACGT') # Illegal VCF can contain '-' to indicate deletions
                    yield SNV.load(self, i)
                else: # 1:N insertion
                    pass
            else:
                if len(alt) == 1: # N:1 deletion
                    pass
                elif len(self.ref) == len(alt): # N:N substitution
                    count = 0
                    unmatched_idx = None
                    for j, (b1, b2) in enumerate(zip(self.ref, alt)):
                        if b1 != b2:
                            count += 1
                            unmatched_idx = j
                    if count == 1: # 1:1 SNV
                        position = self.position + unmatched_idx
                        ref = self.ref[unmatched_idx]
                        alt = self.alts[i][unmatched_idx]
                        AC = self.ACs[i]
                        AF = self.AFs[i]
                        yield SNV(self.chromosome, position, self.rsid, ref, alt, AC, self.AN, AF)
                    else:
                        # Combination of 1:1 SNV (e.g. ACGT->ACCC)
                        raise NotImplementedError
                else: # N:M indel
                    pass

class SNV:
    def __init__(self, chromosome: str, position: int, rsid: str, ref: str, alt: str, AC: int, AN: int, AF: float):
        self.chromosome = normalized_chromosome(chromosome)
        self.position = position
        self.rsid = rsid
        self.ref = ref
        self.alt = alt
        assert(len(ref) == len(alt) == 1)
        assert(ref != alt)

        self.AC = AC
        self.AN = AN
        self.AF = AF

    @classmethod
    def load(cls, entry: VariationEntry, idx: int):
        return cls(entry.chromosome, entry.position, entry.rsid, entry.ref, entry.alts[idx], entry.ACs[idx], entry.AN, entry.AFs[idx])

    def __str__(self) -> str:
        return "{}:{} {}->{} ({}, AF={}={}/{})".format(self.chromosome, self.position, self.ref, self.alt, self.rsid, self.AF, self.AC, self.AN)

def filter_vcf(vcfs, regions: List[Region], fpath, force=False):
    if os.path.exists(fpath) and force is False:
        return

    with open(fpath, 'w') as f:
        for vcf in vcfs:
            args = ["tabix", vcf] # Run external `tabix` command.
            for region in regions:
                args.append("chr{}:{}-{}".format(region.chromosome, region.start, region.end))

            cp = subprocess.run(args, check=True, text=True, capture_output=True)
            if cp.stdout:
                f.write(cp.stdout)

def iterate_vcf(vcf, mode: str) -> Iterator[SNV]:
    assert(mode in ('1KGP', '54KJPN'))

    with open(vcf) as f:
        for l in f:
            try:
                if mode == '1KGP':
                    entry = VariationEntry.load_from_1KGP(l)
                elif mode == '54KJPN':
                    entry = VariationEntry.load_from_54KJPN(l)
                for snv in entry.snvs:
                    yield snv
            except (NotPassedError, BlankLineError):
                continue