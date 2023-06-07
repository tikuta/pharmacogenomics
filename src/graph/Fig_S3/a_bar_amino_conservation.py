#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt

saving_path = '../../../results/graph/Fig_S3/a_bar_amino_conservation'
if not os.path.exists(saving_path):
    os.makedirs(saving_path)
file_error = open('{}/00_error.txt'.format(saving_path),'w')

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
        for g in range(1,len(generic)):
            if len(generic[g]) <= 1:
                continue
            generic_info = list(generic[g])
            amino, pos = generic_info[0], int(''.join(generic_info[1:]))
            GENERIC[gene_name][pos] = GENERIC[gene_name].get(pos, [])
            GENERIC[gene_name][pos].append(ref[g])
            GENERIC[gene_name][pos].append(amino)


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

    MOST[generic] = MOST.get(generic, [])
    MOST[generic].append(max_amino)
    MOST[generic].append(round(max_value/sum_value*100,1))
print(MOST)



def get_score(generic):
    tm = int(generic.split('x')[0])
    num = int(generic.split('x')[1])
    tm_score = tm*10
    if tm > 10:
        tm_score = tm
    num_score = num/100
    if num>100:
        num_score = num/1000
    
    score = tm_score + num_score
    return score

generics = []
values = []
scores = []

for generic in MOST:
    generics.append(generic)
    values.append(MOST[generic][1])
    scores.append(get_score(generic))

generics_ar = np.array(generics)
values_ar = np.array(values)
scores_ar = np.array(scores)

sort = np.argsort(scores)

generics_sort = list(generics_ar[sort])
values_sort = list(values_ar[sort])


fig, ax = plt.subplots(1,1,figsize=(8,2))
keys = ['1x','2x','3x','4x','5x','6x','7x','8x','12x','23x','34x','45x','56x','67x']
colors = ['indigo','blue','royalblue','forestgreen','lime','orange','salmon','red','darkgrey','dimgrey','darkgrey','dimgrey','darkgrey','dimgrey']
for j in range(len(generics_sort)):
    for k in keys:
        if generics_sort[j].startswith(k):
            c = colors[keys.index(k)]
    ax.bar(j,values_sort[j], color=c)

ax.set_xticks(np.arange(len(generics_sort)))
ax.set_xticklabels(generics_sort, rotation=90, fontsize=2)
ax.tick_params(width=0.5, length=2)

fig.savefig('{}/a_bar_amino_conservation.svg'.format(saving_path))