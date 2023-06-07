#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt

saving_path = '../../../results/graph/Fig_3/c_top_AC_variants'
if not os.path.exists(saving_path):
    os.makedirs(saving_path)
file_error = open('{}/00_error.txt'.format(saving_path),'w')


#---------------------------------------------------------------------------------------
#copied from Fig_2/c_bar_structurre_sense.py
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

# Classify the positions according to the structure it belongs to.
ALL_TM_REGION = {}
for gene_name in GENERIC:
    STRUCTURE = {}
    generic_structures = []
    for n in range(1,8):
        key = 'TM{}'.format(str(n))
        STRUCTURE[key] = STRUCTURE.get(key, [])
        generic_structures.append('{}x'.format(str(n)))
    STRUCTURE['H8'] = STRUCTURE.get('H8',[])
    generic_structures.append('8x')
    tms = list(STRUCTURE.keys())


    for pos in GENERIC[gene_name]:
        generic = GENERIC[gene_name][pos][0]
        for i in range(len(generic_structures)):
            if generic.startswith(generic_structures[i]):
                STRUCTURE[tms[i]].append(pos)
    #print(STRUCTURE)
    
    ALL_TM_REGION[gene_name] = ALL_TM_REGION.get(gene_name, {})
    TM_REGION = ALL_TM_REGION[gene_name]
    for tm in tms:
        TM_REGION[tm] = TM_REGION.get(tm, [0,0])
        if len(STRUCTURE[tm]) != 0:
            TM_REGION[tm][0] = min(STRUCTURE[tm])
            TM_REGION[tm][1] = max(STRUCTURE[tm])
    #print(TM_REGION)#inclusive



def which_structure(gene_name, pos):
    TM_REGION = ALL_TM_REGION[gene_name]
    pos = int(pos)
    for tm in tms:
        if TM_REGION[tm][0] <= pos <= TM_REGION[tm][1]:
            structure = tm
            return structure
    if pos < TM_REGION['TM1'][0]:
        structure = 'N'
        return structure
    elif pos > TM_REGION['H8'][1]:
        structure = 'C'
        return structure
    loops = ['ICL1','ECL1','ICL2','ECL2','ICL3','ECL3', 'TM7-8']
    for i in range(8):
        if TM_REGION[tms[i]][1] < pos < TM_REGION[tms[i+1]][0]:
            structure = loops[i]
            return structure
    print(pos)
    return 'error'

structures = ['N','TM1','ICL1','TM2','ECL1','TM3','ICL2','TM4','ECL2','TM5','ICL3','TM6','ECL3','TM7','TM7-8','H8','C']

#---------------------------------------------------------------------------------------

file_isoforms = open('../../../results/3_0_list_isoforms_info/isoforms.txt', 'r')
isoforms = file_isoforms.read().split('\n')
labels = []
numbers = []
for isoform in isoforms:
    if len(isoform) == 0:
        continue
    i = isoform.split('\t')
    is_canonial = i[5]
    if is_canonial != '1':# exclude non canonical isoforms
        continue
    gene_name, m_p, chr_num = i[0], i[3], i[2]
    print(gene_name)
    
    filepath_variant = '../../../results/list/2_0_convert_into_amino/38KJPN/{}.txt'.format(gene_name)
    if not os.path.exists(filepath_variant):
        file_error.write('{}:\tno variant file exists\n'.format(gene_name))
        continue
    print(gene_name)
    file_variant = open(filepath_variant, 'r')
    variants= file_variant.read().split('\n')
    for variant in variants:
        if len(variant) == 0:
            continue
        if not variant.startswith('missense'):
            continue
        v = variant.split('\t')
        if not v[18].startswith('AC'):
            file_error.write('{}_{}:\tinvalid number are written(AC)\n'.format(gene_name, v[18]))
            continue
        if not v[19].startswith('AN'):
            file_error.write('{}_{}:\tinvalid number are written(AF)\n'.format(gene_name, v[19]))
            continue
        pos, before, after, ac, an = v[2], v[7], v[9], int(v[18].split('=')[1]), int(v[19].split('=')[1])
        labels.append('{}_{}{}{}'.format(gene_name, before, pos, after))
        numbers.append(round(ac/an,2))

labels_ar = np.array(labels)
numbers_ar = np.array(numbers)

sort = np.argsort(numbers_ar)

labels_sort = list(labels_ar[sort])
labels_sort.reverse()
numbers_sort = list(numbers_ar[sort])
numbers_sort.reverse()


labels_top = []
numbers_top = []
colors_top = []

structures = ['N','TM1','ICL1','TM2','ECL1','TM3','ICL2','TM4','ECL2','TM5','ICL3','TM6','ECL3','TM7','TM7-8','H8','C']
colors = ['dimgrey','red','darkgrey','salmon','dimgrey','orange','darkgrey','lime','dimgrey','forestgreen','darkgrey','royalblue','dimgrey','blue','darkgrey','indigo','darkgrey']

for i in range(len(numbers_sort)):
    if numbers_sort[i] < 0.5:
        continue
    numbers_top.append(numbers_sort[i])
    labels_top.append(labels_sort[i])
    
    gene_name = labels_sort[i].split('_')[0]
    pos = int(''.join(list(labels_sort[i].split('_')[1])[1:-1]))
    structure = which_structure(gene_name,pos)
    colors_top.append(colors[structures.index(structure)])

labels_top.reverse()
numbers_top.reverse()
colors_top.reverse()

fig, ax = plt.subplots(1,1,figsize=(4,8))
for i in range(len(labels_top)):
    ax.barh(i, numbers_top[i], color = colors_top[i])

ax.set_yticks(np.arange(len(labels_top)))
ax.set_yticklabels(labels_top, fontsize=5)
ax.set_xlabel('Allele Frequency')


labels = ['TM1','TM2','TM3','TM4','TM5','TM6','TM7','H8','Intra Cellular Side','Extra Cellular Side']
colors = ['red','salmon','orange','lime','forestgreen','royalblue','blue','indigo','darkgrey','dimgrey']
for i in range(len(labels)):
    ax.barh(100000, 1000000, color=colors[i], label=labels[i])
ax.set_xlim(0.5,1.03)
ax.set_ylim(-1, len(labels_top)+1)

plt.legend(bbox_to_anchor=(0.4,-0.1), loc='upper center',fontsize=6, ncol=5)

plt.subplots_adjust(left=0.25, right=0.95, bottom=0.3, top=0.95)
fig.savefig('{}/top_AC_missense_variants.svg'.format(saving_path))