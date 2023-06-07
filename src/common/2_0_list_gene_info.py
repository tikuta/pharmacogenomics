#!/usr/bin/env python3

import os
import json

saving_path = '../../results/2_0_list_gene_info'
if not os.path.exists(saving_path):
    os.makedirs(saving_path)
file_error = open('{}/00_error.txt'.format(saving_path),'w')

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
    if gene_name == '':
        file_error.write('{}:\tno gene name\n'.format(line))
        continue
    gpcrs.append(gene_name)
print(gpcrs)


file_write = open('{}/gene.txt'.format(saving_path),'w')
for gpcr in gpcrs:
    print(gpcr)
    filepath_dic = '../../results/1_0_get_json_file/{}.json'.format(gpcr)
    if not os.path.exists(filepath_dic):
        file_error.write('{}:\tno json file\n'.format(gpcr))
        continue
    file_dic = open(filepath_dic, 'r')
    DIC = json.load(file_dic)
    chr_num = DIC['seq_region_name']
    biotype = DIC['biotype']
    emsemble_id = DIC['id']
    strand = DIC['strand']
    start = DIC['start']
    end = DIC['end']
    write = [gpcr, emsemble_id, biotype, chr_num, str(strand), '{}:{}'.format(str(start), str(end))]
    write = '\t'.join(write)
    file_write.write(write + '\n')