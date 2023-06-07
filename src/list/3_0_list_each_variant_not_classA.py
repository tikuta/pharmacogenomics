#!/usr/bin/env python3
#mostly copied from 1_0_list_each_variant_not_classA


import os
import numpy as np

saving_path = '../../results/list/3_0_list_each_variant_not_classA'
if not os.path.exists(saving_path):
    os.makedirs(saving_path)
file_error = open('{}/00_error.txt'.format(saving_path),'w')

file_isoforms = open('../../results/3_0_list_isoforms_info/isoforms.txt', 'r')
isoforms = file_isoforms.read().split('\n')

dataset = '1KGP'
if not os.path.exists('{}/{}'.format(saving_path,dataset)):
    os.makedirs('{}/{}'.format(saving_path,dataset))

for isoform in isoforms:
    if len(isoform) == 0:
        continue
    i = isoform.split('\t')
    is_canonical = i[5]
    if is_canonical != '1':
        continue
    gene_name = i[0]
    
    filepath = '../../results/classify/2_0_class_into_gpcr_not_classA/{}/{}.txt'.format(dataset, gene_name)
    if not os.path.exists(filepath):
        continue
    file_variants = open(filepath, 'r')
    variants = file_variants.read().split('\n')

    file_write = open('{}/{}/{}.txt'.format(saving_path, dataset, gene_name),'w')

    for variant in variants:
        if len(variant) == 0:
            continue

        v = variant.split('\t')
        
        # checked if there are two or more variation at same position in 1KGP
        # concluded that there are no multiple variation in the same low
        # In addition to that, there are no insertion or deletion data in 1KGP dataset
        #if len(v[3]) >= 2:
            #file_error.write('{}\t{}:\tlong variation\n'.format(gene_name, variant))
        write = '\t'.join(v[:8])
        file_write.write(write + '\n')







dataset = '38KJPN'
if not os.path.exists('{}/{}'.format(saving_path,dataset)):
    os.makedirs('{}/{}'.format(saving_path,dataset))

for isoform in isoforms:
    if len(isoform) == 0:
        continue
    i = isoform.split('\t')
    is_canonical = i[5]
    if is_canonical != '1':
        continue
    gene_name = i[0]
    print(gene_name)

    filepath = '../../results/classify/2_0_class_into_gpcr_not_classA/{}/{}.txt'.format(dataset, gene_name)
    if not os.path.exists(filepath):
        continue
    file_variants = open(filepath, 'r')
    variants = file_variants.read().split('\n')

    file_write = open('{}/{}/{}.txt'.format(saving_path, dataset, gene_name),'w')

    for variant in variants:
        if len(variant) == 0:
            continue

        v = variant.split('\t')

        if ',' in v[4]:
            tf = True
            num_variation = len(v[4].split(','))
            write = [[] for _ in range(num_variation)]
        else:
            tf = False
            write = [[]]

        infos = v[7].split(';')

        for info in infos:
            i = info.split('=')
            need = ['AC','AN','AF','nhomalt','nhomalt_XX']
            ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
            ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
            ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency, for each ALT allele, in the same order as listed">
            ##INFO=<ID=nhomalt,Number=A,Type=Integer,Description="Count of homozygous individuals">

            if i[0] in need:
                if tf == False:
                    write[0].append(info)
                    continue
                if i[0] == 'AN':# Since AN doesn't change even if multiple variation exists, only one number(mostly 77444) is written.
                    for k in range(num_variation):
                        write[k].append(info)
                    continue
                i_numbers = i[1].split(',')
                for k in range(num_variation):
                    write[k].append('{}={}'.format(i[0],i_numbers[k]))


        
        for k in range(len(write)):
            final = v[:7] + write[k]
            if tf == True:
                f = final[4].split(',')
                final[4] = f[k]
            final = '\t'.join(final)
            file_write.write(final + '\n')
