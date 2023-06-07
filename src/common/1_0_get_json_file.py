#!/usr/bin/env python3

import os
import requests
import json

saving_path = '../../results/1_0_get_json_file'
if not os.path.exists(saving_path):
    os.makedirs(saving_path)
file_error = open('{}/00_error.txt'.format(saving_path), 'w')

file_gpcr = open('../../data/NCBi/GPCRTargets.csv','r')
lines = file_gpcr.read().split('\n')
gpcrs = []
for line in lines:
    if len(line) <= 1:
        continue
    if line.startswith('"Type"'):
        continue
    l = line.split('","')
    gene_name = l[10]
    gpcrs.append(gene_name)


for gpcr in gpcrs:
    print(gpcr)
    filepath_write = '{}/{}.json'.format(saving_path,gpcr)
    if os.path.exists(filepath_write):
        continue
    file_write = open(filepath_write, 'w')

    server = server = "https://rest.ensembl.org"
    ext = "/lookup/symbol/homo_sapiens/{}?expand=1".format(gpcr)
    # To explore on internet, search for 'https://rest.ensembl.org/lookup/symbol/homo_sapiens/{}?content-type=application/json;expand=1'.format(gpcr)
    # The sorce script is on this website 'https://rest.ensembl.org/documentation/info/symbol_lookup'
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    
    if not r.ok:
        res = r.status_code
        print(res)
        file_error.write(gpcr + ':\t' + str(res) + '\n')
        continue

    decoded = r.json()
    #print(repr(decoded))

    json.dump(decoded, file_write, ensure_ascii=False)
    