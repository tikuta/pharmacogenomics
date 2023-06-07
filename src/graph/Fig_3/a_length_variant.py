#!/usr/bin/env python3

import os
import matplotlib.pyplot as plt

saving_path = '../../../results/graph/Fig_3/a_length_variant'
if not os.path.exists(saving_path):
    os.makedirs(saving_path)
file_error = open('{}/00_error.txt'.format(saving_path),'w')


#make dictionary (key = common name, value = gene name)
file_name = open('../../../data/GPCRdb/handmade_gene_name.csv','r')
names = file_name.read().split('\n')
NAME = {}
for name in names:
    n = name.split(',')
    NAME[n[0].replace('\ufeff','')] = NAME.get(n[0].replace('\ufeff',''),'')
    NAME[n[0].replace('\ufeff','')] = n[1] 
# add gpcrs belonging to other classes to the 'NAME' dictionary
file_name = open('../../../data/GPCRdb/handmade_gene_name_other_classes.csv','r')
names = file_name.read().split('\n')
for name in names:
    n = name.split(',')
    NAME[n[0].replace('\ufeff','')] = NAME.get(n[0].replace('\ufeff',''),'')
    NAME[n[0].replace('\ufeff','')] = n[1] 

# Make dictionary( key = gene name, value = class)
gpcr_classes = ['A','B1_secretin', 'B2_Adhesion','C_Glutamate','F_Frizzled','T_Taste']
CLASS = {}
for gpcr_class in gpcr_classes:
    file_alignment = open('../../../data/GPCRdb/GPCRdb_alignment_{}.csv'.format(gpcr_class),'r')
    lines = file_alignment.read().split('\n')
    for line in lines:
        if len(line) == 0:
            continue
        if not line.startswith('[Human]'):
            continue
        receptor = line.split(',')[0].replace('[Human]','').replace('receptor','').replace(' ','')
        if receptor in NAME:
            receptor = NAME[receptor]
        CLASS[receptor] = CLASS.get(receptor,'')
        CLASS[receptor] = gpcr_class
print(CLASS)






def give_path(gene_name):
    if not gene_name in CLASS:
        return 'error'
    if CLASS[gene_name] == 'A':
        return '../../../results/list/2_0_convert_into_amino/38KJPN/{}.txt'.format(gene_name)
    return '../../../results/list/4_0_convert_into_amino_not_classA/38KJPN/{}.txt'.format(gene_name)


def get_length(positions):
    count = 0
    for position in positions:
        pos = position.split(':')
        count = count + (int(pos[1]) - int(pos[0]) + 1)
    amino_count = count // 3
    return count, amino_count



senses = ['missense','silent','nonsense']
count_all = []

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
    positions = i[6:]
    print(gene_name)
    length = list(get_length(positions))
    if length[0] % 3 != 0:
        file_error.write('{}:\tincorrect base length\n'.format(gene_name))
    print(length)

    filepath_variant = give_path(gene_name)
    if filepath_variant == 'error':
        file_error.write('{}:\tno variant information\n'.format(gene_name))
        continue
    file_variant = open(filepath_variant,'r')
    variants = file_variant.read().split('\n')
    count = [0 for _ in range(len(senses))]
    for variant in variants:
        if len(variant) == 0:
            continue
        sense = variant.split('\t')[0]

        if sense == 'insertion' or sense == 'deletion':
            continue
        
        count[senses.index(sense)] = count[senses.index(sense)] + 1
    
    #only consider the missense variant(count[0])
    count_all.append([gene_name, length[0], length[1], count[0]])
print(count_all)


fig = plt.figure()

ax1 = fig.subplots(1,1)
ax2 = ax1.twiny()

colors = ['tomato','darkgreen','limegreen','royalblue','aqua','purple']
markers = ['o','^','^','^','^','^']
file_statistic = open('{}/raw_data.csv'.format(saving_path),'w')
file_statistic.write('gene name\tbase length\tamino acid length\tnumber of missense variant\n')
for i in range(len(count_all)):
    n = gpcr_classes.index(CLASS[count_all[i][0]])
    ax1.scatter([count_all[i][2]], [count_all[i][3]], color=colors[n], marker=markers[n], s=5)
    ax2.scatter([count_all[i][1]], [count_all[i][3]], alpha=0)

ax1.set_xlabel('Protein Length (residues)')
ax2.set_xlabel('Protein Length (bases)')
ax1.set_ylabel('Number of Missense Variants')

labels = ['A','B1(secretin)','B2(Adhesion)','C','F','Taste']
for j in range(len(labels)):
    ax1.scatter(-1000, -1000, color = colors[j], label = labels[j], marker = markers[j])
ax1.set_xlim(0,)
ax1.set_ylim(0,)
ax1.legend(bbox_to_anchor=(0.5,-0.2), loc='upper center',fontsize=6, ncol=3)
plt.subplots_adjust(left=0.2, right=0.95, bottom=0.3, top=0.9)
fig.savefig('{}/a_length_correlation_missense.svg'.format(saving_path))