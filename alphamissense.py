#!/usr/bin/env python3
import gzip
import os
from typing import List
import gpcrdb
from config import AM_FILENAME
from vcf import SNV
from utils import normalized_chromosome

def _extract(accessions: List[str], paths: List[str]):
    assert(len(accessions) == len(paths))

    am = gzip.open(os.path.join("data", "AlphaMissense_hg38.tsv.gz"), 'rt')
    for l in am:
        if l.startswith('#CHROM'): # header
            for path in paths:
                with open(path, 'w') as f:
                    f.write(l)
        elif l.startswith('#'):
            continue
        else:
            cols = l.strip().split('\t')
            accession = cols[5]
            if accession in accessions:
                print(cols[0], cols[1], end='\r')
                with open(paths[accessions.index(accession)], 'a') as f:
                    f.write(l)

def extract(force=False):
    accessions, paths = [], []
    for gpcrdb_entry in gpcrdb.get_filtered_receptor_list(force=force):
        accession = gpcrdb_entry.accession
        path = os.path.join(gpcrdb_entry.dirpath, AM_FILENAME)

        if force is True or not os.path.exists(path):
            accessions.append(accession)
            paths.append(path)
    print("Total {} genes to be extracted from AlphaMissense dataset.".format(len(accessions)))
    
    if len(accessions) > 0:
        _extract(accessions, paths)

def _lookup(path, snv: SNV) -> float:
    with open(path) as f:
        for l in f:
            if l.startswith('#'):
                continue
            cols = l.strip().split('\t')
            position = int(cols[1])
            ref, alt = cols[2], cols[3]
            uniprot_id, transcript_id = cols[5], cols[6]
            protein_variant = cols[7]
            am_pathogenicity, am_class = float(cols[8]), cols[9]
            
            if snv.position == position:
                assert(snv.ref == ref)
                if snv.alt == alt:
                    return am_pathogenicity
    raise Exception("No variation found in AlphaMissense ({}).".format(snv))

def _check_before_lookup(cols: List[str], expected_chromosome: str, expected_uniprot_id: str, expected_transcript_id: str):
    chromosome = normalized_chromosome(cols[0])
    genome = cols[4]
    uniprot_id = cols[5]
    transcript_id = cols[6]
    assert(genome == 'hg38')
    assert(chromosome == expected_chromosome)
    assert(uniprot_id == expected_uniprot_id)
    assert(transcript_id == expected_transcript_id)

def _check_after_lookup(cols: List[str], expected_substitution: str):
    substitution = cols[7]
    assert(substitution == expected_substitution)

def lookup(path, snv: SNV, expected_uniprot_id: str, expected_transcript_id: str, expected_substitution: str):
    with open(path) as f:
        for ln, l in enumerate(f):
            if l.startswith('#'):
                continue

            cols = l.strip().split('\t')
            if ln == 1:
                _check_before_lookup(cols, snv.chromosome, expected_uniprot_id, expected_transcript_id)
            
            position = int(cols[1])
            ref, alt = cols[2], cols[3]
            am_pathogenicity = float(cols[8])
            
            if snv.position == position and snv.ref == ref and snv.alt == alt:
                _check_after_lookup(cols, expected_substitution)
                return am_pathogenicity
    raise Exception("No variation found in AlphaMissense ({}).".format(snv))