#!/usr/bin/env python3

import os
import requests
import json

save_dir = './lookup'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

file_gpcr = open('../GtoP/GPCRTargets.csv','r')
lines = file_gpcr.read().split('\n')
gpcrs = set() # there are some duplicates in the list (e.g. GPR42)
for line in lines:
    if len(line) <= 1:
        continue
    if line.startswith('"Type"'):
        continue
    l = line.split('","')
    gene_name = l[10] # HGNC symbol
    if len(gene_name) > 0: # some genes are not present in human (only in rodents)
        gpcrs.add(gene_name)

for gpcr in gpcrs:
    filepath_write = os.path.join(save_dir, '{}.json'.format(gpcr))
    file_write = open(filepath_write, 'w')

    server = "https://rest.ensembl.org"
    ext = "/lookup/symbol/homo_sapiens/{}?expand=1".format(gpcr)
    uri = server + ext
    print(uri)

    r = requests.get(uri, headers={"Content-Type": "application/json"})

    if not r.ok:
        raise Exception

    json.dump(r.json(), file_write, ensure_ascii=False)
