#!/usr/bin/env python3
# This analysis wasn't needed

import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

saving_path = '../../../results/graph/Fig_3/b_1_variant_AC_accumulation'
if not os.path.exists(saving_path):
    os.makedirs(saving_path)
file_error = open('{}/00_error.txt'.format(saving_path),'w')

#---------------------------------------------------------------------------------------------------------------

#make dictionary (key = common name, value = gene name)
file_name = open('../../../data/GPCRdb/handmade_gene_name_2.csv','r')
names = file_name.read().split('\n')
NAME = {}
for name in names:
    n = name.split(',')
    NAME[n[0].replace('\ufeff','')] = NAME.get(n[0].replace('\ufeff',''),'')
    NAME[n[0].replace('\ufeff','')] = n[1] 

#get generic number
GENERIC = {}
GENERIC_LS = {}
filepath_generics = ['../../../data/GPCRdb/generic_number/residue_table_startswithO.csv','../../../data/GPCRdb/generic_number/residue_table_startswithP.csv','../../../data/GPCRdb/generic_number/residue_table_startswithQ.csv']
for filepath_generic in filepath_generics:
    file_generic = open(filepath_generic)
    contents = file_generic.read().split('\n')
    generic_list = []
    for n in range(len(contents)):
        content = contents[n]
        c = content.split(',')
        generic_list.append(c)
    generic_ar = np.array(generic_list)
    generic_list = np.transpose(generic_ar)
    for generic in generic_list:
        if generic[0] == '\ufeffGPCRdb(A)':
            ref = generic
            continue
        if generic[0] == 'GPCRdb(A)' or generic[0] == 'BW':
            continue
        gene_name = generic[0].split(' ')[0].replace('<i>','').replace('</i>','')
        if gene_name in NAME:
            gene_name = NAME[gene_name]
        GENERIC[gene_name] = GENERIC.get(gene_name, {})
        GENERIC_LS[gene_name] = GENERIC_LS.get(gene_name, [])
        for g in range(1,len(generic)):
            if len(generic[g]) <= 1:
                continue
            generic_info = list(generic[g])
            amino, pos = generic_info[0], int(''.join(generic_info[1:]))
            GENERIC[gene_name][pos] = GENERIC[gene_name].get(pos, [])
            GENERIC[gene_name][pos].append(ref[g])
            GENERIC[gene_name][pos].append(amino)
            GENERIC_LS[gene_name].append(ref[g])

print(GENERIC_LS.keys())


CONSERVE = {}
for gene_name in GENERIC:
    #exclude not class A GPCR
    filepath = '../../../results/list/2_0_convert_into_amino/38KJPN/{}.txt'.format(gene_name)
    if not os.path.exists(filepath):
        continue
    for pos in GENERIC[gene_name]:
        generic = GENERIC[gene_name][pos][0]
        amino = GENERIC[gene_name][pos][1]
        
        CONSERVE[generic] = CONSERVE.get(generic, {})
        CONSERVE[generic][amino] = CONSERVE[generic].get(amino,0)
        CONSERVE[generic][amino] = CONSERVE[generic][amino] + 1


MOST = {}
for generic in CONSERVE:
    DIC = CONSERVE[generic]
    aminos = list(DIC.keys())
    values = list(DIC.values())
    sum_value = sum(values)
    if sum_value < len(list(GENERIC.keys()))*0.5:#only consider the position that is named generic number in more than 50% of class A GPCRs.
        continue
    '''
    print(generic)
    print(sum_value)
    print(aminos)
    print(values)
    '''
    max_value = max(values)
    max_amino = aminos[values.index(max_value)]

    MOST[generic] = MOST.get(generic, 0)
    MOST[generic] = sum_value
#print(MOST)
conserved_generics = list(MOST.keys())
#print(conserved_generics)

#-----------------------------------------------------------------------------------------------------------------------
def get_score(generic):
    tm = int(generic.split('x')[0])
    num = int(generic.split('x')[1])
    tm_score = tm*10
    if tm > 10:
        tm_score = tm
    num_score = num/100
    if num > 100:
        num_score = num/1000
    
    score = tm_score + num_score
    return score


scores = []
for generic in conserved_generics:
    score = get_score(generic)
    scores.append(score)

scores_ar = np.array(scores)
conserved_generics_ar = np.array(conserved_generics)

sort = np.argsort(scores_ar)

conserved_generics = list(conserved_generics_ar[sort])
#print(conserved_generics)


# Count the number of varants that exist on the specific position-------------------------------------------------------
COUNT = {}
count_gpcrs = 0

file_isoforms = open('../../../results/3_0_list_isoforms_info/isoforms.txt', 'r')
isoforms = file_isoforms.read().split('\n')
for isoform in isoforms:
    if len(isoform) == 0:
        continue
    i = isoform.split('\t')
    is_canonial = i[5]
    if is_canonial != '1':# exclude non canonical isoforms
        continue
    gene_name, m_p, chr_num = i[0], i[3], i[2]



    filepath_variat = '../../../results/list/2_0_convert_into_amino/38KJPN/{}.txt'.format(gene_name)
    if not os.path.exists(filepath_variat):
        file_error.write('{}:\tno variant file exists\n'.format(gene_name))
        continue
    count_gpcrs = count_gpcrs + 1
    file_variat = open(filepath_variat,'r')
    variants = file_variat.read().split('\n')



    COUNT[gene_name] = COUNT.get(gene_name, [0 for _ in range(len(conserved_generics))])
    if not gene_name in GENERIC_LS:
        file_error.write('{}:\tnot proper gene_name are in GENERIC_LS or the gene_name not exists in the dictionary\n'.format(gene_name))
        continue
    isoform_generics = GENERIC_LS[gene_name]
    for  generic in conserved_generics:
        if not generic in isoform_generics:
            index = conserved_generics.index(generic)
            COUNT[gene_name][index] = 0


    for variant in variants:
        if not variant.startswith('missense'):
            continue
        v = variant.split('\t')
        generic = v[5]
        if generic == '.':
            continue
        if not generic in conserved_generics:
            continue
        index = conserved_generics.index(generic)
        info = v[17:]
        for i in info:
            if i.startswith('AF='):
                af = float(i.split('=')[1])
                COUNT[gene_name][index] = af

#print(COUNT)
#print(count_gpcrs)





#------------------------------------------------------------------------------------------------------------------------


gpcrs = np.array(list(COUNT.keys()))
gpcrs = np.sort(gpcrs)

values = []
for gpcr in gpcrs:
    values.append(COUNT[gpcr])

values = np.array(values)
print(values)


fig, ax = plt.subplots(1,1)
sns.heatmap(data=values, cmap='Reds')
fig.savefig('{}/heatmap2.png'.format(saving_path))