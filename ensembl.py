#!/usr/bin/env python3
from typing import List
import requests
import json
import os
from utils import Region

def get_gene(gene_id, fpath, force=False):
    if os.path.exists(fpath) and force is False:
        with open(fpath) as f:
            return json.load(f)
        
    uri = "https://rest.ensembl.org/lookup/id/{}?expand=1".format(gene_id)
    r = requests.get(uri, headers={"Content-Type": "application/json"})

    if not r.ok:
        raise Exception
    
    j = r.json()
    
    with open(fpath, 'w') as f:
        json.dump(j, f, indent=2)

    return j


def get_sequence(region: Region, strand: int, fpath, force=False):
    assert(strand == 1 or strand == -1)
    if os.path.exists(fpath) and force is False:
        with open(fpath) as f:
            return json.load(f)
        
    uri = "https://rest.ensembl.org/sequence/region/human"

    data = {
        "coord_system_version": "GRCh38",
        "regions": ["{}:{}..{}:{}".format(region.chromosome, region.start, region.end, strand)]
    }
    r = requests.post(uri, headers={"Content-Type": "application/json"}, data=json.dumps(data))

    if not r.ok:
        r.raise_for_status()
        raise Exception
    
    j = r.json()

    with open(fpath, 'w') as f:
        json.dump(j, f, indent=2)
        
    return j

def get_coding_regions(transcript) -> List[Region]:
    assert('Translation' in transcript)

    chromosome = transcript['seq_region_name']

    # We only have positions where coding regions start and end
    #                    and where exons          start and end.
    # So we must assembly exact coding regions (indicated as '*' in the following examples).
    # Note that start (s) is smaller than end (e), reagrdless of strand (plus or minus).
    # 
    #                   cds_s                                cds_e
    #                     |                                    |
    # ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
    #                     *****       *********       **********
    #  |---|     |------------|       |-------|       |--------------|        |---| 
    #  exon1   exon2_s      exon2_e exon3_s exon3_e exon4_s        exon4_e    exon5
    #
    # Or,
    #                   cds_s                                 cds_e
    #                     |                                     |
    # ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
    #                     ***************************************
    #             |------------------------------------------------------------|
    #           exon6_s                                                      exon6_e

    cds_s = transcript['Translation']['start']
    cds_e = transcript['Translation']['end']
    assert(cds_s < cds_e)

    coding_regions = []
    for exon in transcript['Exon']:
        exon_s = int(exon['start'])
        exon_e = int(exon['end'])
        assert(exon_s < exon_e)
        
        if exon_s <= cds_s <= exon_e <= cds_e:   # Exon 2 in the above example
            coding_regions.append(Region(chromosome, cds_s, exon_e))
        elif cds_s <= exon_s <= exon_e <= cds_e: # Exon 3
            coding_regions.append(Region(chromosome, exon_s, exon_e))
        elif cds_s <= exon_s <= cds_e <= exon_e: # Exon 4
            coding_regions.append(Region(chromosome, exon_s, cds_e))
        elif exon_s <= cds_s <= cds_e <= exon_e: # Exon 6
            coding_regions.append(Region(chromosome, cds_s, cds_e))
        else:                                    # Exons 1 and 5 contain no coding region
            continue

    return coding_regions


