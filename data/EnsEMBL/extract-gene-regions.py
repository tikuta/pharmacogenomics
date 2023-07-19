#!/usr/bin/env python3

import os
import json
import glob

with open("gene-regions.txt", 'w') as f:
    for p in glob.glob("./lookup/*.json"):
        gene = os.path.splitext(os.path.basename(p))[0]
        j = json.load(open(p, 'r'))

        biotype = j['biotype']
        if biotype != "protein_coding":
            continue

        chr_num = j['seq_region_name']
        ensg_id = j['id']
        strand = j['strand'] # 1 or -1
        start = j['start']
        end = j['end']
        write = [gene, ensg_id, "chr" + chr_num, str(strand), '{}-{}'.format(start, end)]
        write = '\t'.join(write)
        f.write(write + '\n')