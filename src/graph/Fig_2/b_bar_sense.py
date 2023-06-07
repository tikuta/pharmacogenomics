#!/usr/bin/env python3

import os
import numpy as np
import sys
import matplotlib.pyplot as plt

saving_path = '../../../results/graph/Fig_2/b_bar_sense'
if not os.path.exists(saving_path):
    os.makedirs(saving_path)
file_error = open('{}/00_error.txt'.format(saving_path),'w')


def count(line, c):
    l = line.split('\t')
    before, after = l[3], l[4]

    snv = False
    indel = False

    if len(before) == len(after) and len(after) == 1:
        c[0] = c[0] + 1
        c[1] = c[1] + 1

    elif ',' in after:
        a_ls = after.split(',')
        c[0] = c[0] + 1
        for a in a_ls:
            if len(a) == len(before):
                if snv == False:
                    c[1] = c[1] + 1
                    snv = True
                
            else:
                if indel == False:
                    c[2] = c[2] + 1
                    indel = True
    
    elif len(before) != len(after):
        c[0] = c[0] + 1
        c[2] = c[2] + 1


senses = ['missense','silent','nonsense','insertion','deletion']
count_sense = [0,0,0,0,0]

file_stats = open('{}/row_data.csv'.format(saving_path),'w')
file_stats.write('gene_name\t'+'\t'.join(senses) + '\n')

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

    filepath_variant = '../../../results/list/2_0_convert_into_amino/38KJPN/{}.txt'.format(gene_name)
    if not os.path.exists(filepath_variant):
        file_error.write('{}:\tvariant file path does not exist\n')
        continue
    file_variant = open(filepath_variant,'r')
    variants = file_variant.read().split('\n')
    count_sense_each = [0,0,0,0,0]
    for variant in variants:
        if len(variant) == 0:
            continue
        info = variant.split('\t')
        sense = info[0]
        index = senses.index(sense)
        count_sense_each[index] = count_sense_each[index] + 1
    
    str_count_sense_each = []
    for c in range(len(count_sense_each)):
        str_count_sense_each.append(str(count_sense_each[c]))
        count_sense[c] = count_sense[c] + count_sense_each[c]
    
    file_stats.write(gene_name + '\t' + '\t'.join(str_count_sense_each) + '\n')

str_count_sense = []
for c in count_sense:
    str_count_sense.append(str(c))
file_stats.write('total\t'+'\t'.join(str_count_sense)+ '\n')

#indel included
#percent_count_sense = [round(count_sense[0]/sum(count_sense)*100, 1), round(count_sense[1]/sum(count_sense)*100, 1),round(sum(count_sense[2:])/sum(count_sense)*100, 1)]

#indel excluded
percent_count_sense = [round(count_sense[0]/sum(count_sense[:3])*100, 1), round(count_sense[1]/sum(count_sense[:3])*100, 1),round(count_sense[2]/sum(count_sense[:3])*100, 1)]

fig, ax = plt.subplots(1,1, figsize=(4,1))
colors = ['lightgrey','grey','tomato']
for i in range(len(percent_count_sense)):
    ax.barh(0, percent_count_sense[i], left=sum(percent_count_sense[:i]), color=colors[i])

labels = ['missense\n{}%'.format(str(percent_count_sense[0])), 'silent\n{}%'.format(str(percent_count_sense[1])), 'nonsense\n{}%'.format(str(percent_count_sense[2]))]
ax.text(percent_count_sense[0]/2, 0, labels[0], ha='center', va='center', fontfamily='Arial')
for i in range(1,len(percent_count_sense)):
    ax.text(percent_count_sense[i]/2 + sum(percent_count_sense[:i]), 0, labels[i], ha='center', va='center', fontfamily='Arial')

#tickを消す
ax.set_xticks(np.arange(0))
ax.set_yticks(np.arange(0))
#囲い線を消す
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)



fig.savefig('{}/sense.svg'.format(saving_path))