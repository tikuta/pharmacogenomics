#!/usr/bin/env python3
# This src is to clarify the reason why GPR42 exists more than one in the lst
# It was because of the GuideToPharmacology's list (/data/NCBi/GPCRTargets.csv)

import os 

saving_path = '../../results/list/2_1_check'
if not os.path.exists(saving_path):
    os.makedirs(saving_path)
file_error = open('{}/00_error.txt'.format(saving_path), 'w')

file_result = open('{}/result.txt'.format(saving_path), 'w')

file_isoforms = open('../../results/3_0_list_isoforms_info/isoforms.txt', 'r')
isoforms = file_isoforms.read().split('\n')
gene_names = []
for isoform in isoforms:
    if len(isoform) == 0:
        continue
    i = isoform.split('\t')
    is_canonial = i[5]
    if is_canonial != '1':# exclude non canonical isoforms
        continue
    gene_name, m_p, chr_num = i[0], i[3], i[2]
    if gene_name in gene_names:
        file_result.write('{}\texist more than two\n'.format(gene_name))
    gene_names.append(gene_name)

    filepath_variant = '../../results/list/2_0_convert_into_amino/38KJPN/{}.txt'.format(gene_name)
    if not os.path.exists(filepath_variant):
        file_error.write('{}:\tno variant file exists\n'.format(gene_name))
        continue

    file_variant = open(filepath_variant,'r')
    variants = file_variant.read().split('\n')

    variant_ls = []
    tf = False
    for variant in variants:
        if not variant in variant_ls:
            variant_ls.append(variant)
            continue
        tf = True
    if tf == True:
        file_result.write(gene_name + '\n')