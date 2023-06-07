#!/usr/bin/env python3

import os
import json
import glob

saving_path = '../../results/2_gene-regions'
if not os.path.exists(saving_path):
    os.makedirs(saving_path)

file_write = open('{}/genes.txt'.format(saving_path),'w')
for p in glob.glob("../../results/1_genes/*.json"):
    gene = os.path.splitext(os.path.basename(p))[0]
    j = json.load(open(p, 'r'))
    chr_num = j['seq_region_name']
    biotype = j['biotype']
    emsemble_id = j['id']
    strand = j['strand'] # 1 or -1
    start = j['start']
    end = j['end']
    write = [gene, emsemble_id, biotype, chr_num, str(strand), '{}:{}'.format(start, end)]
    write = '\t'.join(write)
    file_write.write(write + '\n')