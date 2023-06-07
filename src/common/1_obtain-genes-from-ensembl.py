#!/usr/bin/env python3

import os
import requests
import json

saving_path = '../../results/1_genes'
if not os.path.exists(saving_path):
    os.makedirs(saving_path)

file_gpcr = open('../../data/GtoP/GPCRTargets.csv','r')
lines = file_gpcr.read().split('\n')
gpcrs = set() # there are some duplicates in the list (e.g. GPR42)
for line in lines:
    if len(line) <= 1:
        continue
    if line.startswith('"Type"'):
        continue
    l = line.split('","')
    gene_name = l[10] # HGNC symbol
    if len(gene_name) > 0: # some genes are not present in the human genome (only in rats)
        gpcrs.add(gene_name)

for gpcr in gpcrs:
    print(gpcr)
    filepath_write = '{}/{}.json'.format(saving_path, gpcr)
    file_write = open(filepath_write, 'w')

    server = "https://rest.ensembl.org"
    ext = "/lookup/symbol/homo_sapiens/{}?content-type=application/json;expand=1".format(gpcr)

    r = requests.get(server + ext)
    
    if not r.ok:
        raise Exception

    json.dump(r.json(), file_write, ensure_ascii=False)
