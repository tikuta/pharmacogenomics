#!/usr/bin/env python3

import os
import gzip

saving_path = '../../results/filter/1_0_fileter_vcf_on_GPCR_gene'
if not os.path.exists(saving_path):
    os.makedirs(saving_path)
file_error = open('{}/00_error.txt'.format(saving_path),'w')

#this gene_list will tell you the informatio of chromosome number, start position of gpcr gene, and end position of gpcr gene.
file_gene = open('../../results/2_0_list_gene_info/gene.txt','r')
genes = file_gene.read().split('\n')
genes_list = []
for gene in genes:
    #print(gene)
    if len(gene) == 0 :
        continue
    chr_m = gene.split('\t')[3]
    gs = gene.split('\t')[5].split(':')
    #print(gs)
    region = [chr_m]
    for g in gs:
        region.append(int(g))
    genes_list.append(region)

#print(genes_list)

def check_on_gpcr(c_num, pos):
    tf = False
    for genes in genes_list:
        if c_num != genes[0]:
            continue
        if genes[1] <= pos <= genes[2]:
            tf = True
            break
        else:
            tf = False
    if tf == True:
        return True
    else:
        return False


datasets = ['1KGP','38KJPN']

for dataset in datasets:
    if dataset == '38KJPN':
        filepaths = ['../../data/38KJPN/tommo-38kjpn-20220929-GRCh38-gf-chrX_PAR2.vcf.gz','../../data/38KJPN/tommo-38kjpn-20220929-GRCh38-gf-chrX_PAR3.vcf.gz','../../data/38KJPN/tommo-38kjpn-20220630-GRCh38-gf-autosome.vcf.gz']
        continue
    if dataset == '1KGP':
        filepaths = ['../../data/1KGP/ALL.chrX.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz']
        for c in range(1,23):
            filepaths.append('../../data/1KGP/ALL.chr{}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz'.format(str(c)))
    #filepaths = ['../../data/1KGP/ALL.chr8.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz']
    saving_path_dataset = '{}/{}'.format(saving_path,dataset)
    if not os.path.exists(saving_path_dataset):
        os.makedirs(saving_path_dataset)
    file_info = open('{}/info.txt'.format(saving_path_dataset),'w')
    file_write = open('{}/on_GPCR_gene.vcf'.format(saving_path_dataset),'w')
    for filepath in filepaths:
        print(filepath)
        file = gzip.open(filepath, 'rt')
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
            
            if count % 100000 == 0:
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