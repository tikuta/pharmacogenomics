#!/usr/bin/env python3

import os
import gzip

saving_path = '../../results/filter/2_0_fileter_vcf_on_canonical_isoform'
if not os.path.exists(saving_path):
    os.makedirs(saving_path)
file_error = open('{}/00_error.txt'.format(saving_path),'w')

file_gene = open('../../results/3_0_list_isoforms_info/isoforms.txt','r')
genes = file_gene.read().split('\n')
all_region_list = []
for gene in genes:
    if len(gene) == 0 :
        continue
    is_canonical = gene.split('\t')[5]
    if is_canonical != '1':
        continue
    chr_m = gene.split('\t')[2]
    regions = gene.split('\t')[6:]
    #print(gs)
    region_list = [chr_m]
    for region in regions:
        rs = region.split(':')
        ls = []
        for r in rs:
            ls.append(int(r))
        region_list.append(ls)
    all_region_list.append(region_list)


def check_on_gpcr(c_num, pos):
    for region_list in all_region_list:
        if c_num != region_list[0]:
            tf = False
            continue
        regions = region_list[1:]
        #print(regions)
        tf = False
        for region in regions:
            #print(region)
            if region[0] <= pos <= region[1]:
                tf = True
                break
            else:
                tf = False
        if tf == True:
            break
    if tf == True:
        return True
    else:
        return False


datasets = ['1KGP','38KJPN']

for dataset in datasets:
    saving_path_dataset = '{}/{}'.format(saving_path,dataset)
    if not os.path.exists(saving_path_dataset):
        os.makedirs(saving_path_dataset)
    file_info = open('{}/info.txt'.format(saving_path_dataset),'w')
    file_write = open('{}/on_canonical_isoforms.vcf'.format(saving_path_dataset),'w')
    file_path = '../../results/filter/1_0_fileter_vcf_on_GPCR_gene/{}/on_GPCR_gene.vcf'.format(dataset)
    file = open(file_path, 'r')
    line = file.readline()
    count = 1
    while line:
        if line.startswith('#'):
            file_info.write(line)
            line = file.readline()
            continue
        if len(line) == 0:
            line = file.readline()
            continue
        
        if count % 10000 == 0:
            print(count)
        count = count + 1
        
        l = line.split(('\t'))
        c_num = l[0].replace('chr','')
        pos = int(l[1])
        on_gpcr = check_on_gpcr(c_num,pos)
        if on_gpcr == False:
            line = file.readline()
            continue
        file_write.write(line)
        line = file.readline()