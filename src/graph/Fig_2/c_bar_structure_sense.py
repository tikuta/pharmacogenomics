#!/usr/bin/env python3

import os
import matplotlib.pyplot as plt
import numpy as np

saving_path = '../../../results/graph/Fig_2/c_bar_structure_sense'
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



def which_structure(pos):
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



# Count the total number of variant in each structure.
file_statistic = open('{}/raw_data.csv'.format(saving_path),'w')

structures = ['N','TM1','ICL1','TM2','ECL1','TM3','ICL2','TM4','ECL2','TM5','ICL3','TM6','ECL3','TM7','TM7-8','H8','C']
senses = ['missense','silent','nonsense']
file_statistic.write('.\t{}\n'.format('\t'.join(structures)))
count_all = [[0 for _ in range(len(structures))] for _ in range(3)]
print(count_all)

file_isoforms = open('../../../results/3_0_list_isoforms_info/isoforms.txt', 'r')
isoforms = file_isoforms.read().split('\n')
for isoform in isoforms:
    if len(isoform) == 0:
        continue
    i = isoform.split('\t')
    is_canonial = i[5]
    if is_canonial != '1':# exclude non canonical isoforms
        continue
    gene_name = i[0]

    print(gene_name)

    count_each = [[0 for _ in range(len(structures))] for _ in range(3)]

    filepath_variant = '../../../results/list/2_0_convert_into_amino/38KJPN/{}.txt'.format(gene_name)
    if not os.path.exists(filepath_variant):
        file_error.write('{}:\t no variant file\n'.format(gene_name))
        continue
    file_variant = open(filepath_variant,'r')
    variants = file_variant.read().split('\n')
    for variant in variants:
        if len(variant) == 0:
            continue
        sense, pos = variant.split('\t')[0],variant.split('\t')[2]
        if sense == 'insertion' or sense == 'deletion':
            continue
        structure_name = which_structure(pos)
        index = structures.index(structure_name)
        sense_index = senses.index(sense)
        count_all[sense_index][index] = count_all[sense_index][index] + 1
        count_each[sense_index][index] = count_each[sense_index][index] + 1
    
    file_statistic.write('{}\t'.format(gene_name))
    for j in range(len(senses)):
        for k in range(len(structures)):
            file_statistic.write(str(count_each[j][k]))
            file_statistic.write('\t')
    file_statistic.write('\n')
print(count_all)

file_statistic.write('total\t')
for j in range(len(senses)):
    for k in range(len(structures)):
        file_statistic.write(str(count_all[j][k]))
        file_statistic.write('\t')
file_statistic.write('\n')


fig, ax = plt.subplots(1,1, figsize=(4,4))

width=0.3

xs = np.arange(len(structures))
colors = ['tomato','royalblue','limegreen']
for j in range(len(senses)):
    for k in range(len(structures)):
        position = xs[k] - width + width*j
        ax.bar(position, count_all[j][k], width=width, color=colors[j])


labels = ['missense','silent','nonsense']
for j in range(len(labels)):
    ax.bar(1000, 1000, color=colors[j], label=labels[j])

ax.set_xticks(np.arange(len(structures)))
ax.set_xticklabels(structures, rotation=90)
ax.set_xlim(-1,len(structures))

plt.legend(bbox_to_anchor=(0.5,-0.2), loc='upper center',fontsize=6, ncol=3)
plt.subplots_adjust(left=0.25, right=0.95, bottom=0.3, top=0.95)
fig.savefig('{}/bar_structure_sense.svg'.format(saving_path))