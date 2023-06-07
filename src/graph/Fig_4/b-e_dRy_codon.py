#!/usr/bin/env python3
# 図を作成したものの、示し方を検討中


import os
import numpy as np
import matplotlib.pyplot as plt

saving_path = '../../../results/graph/Fig_4/b_e_dRy_codon'
if not os.path.exists(saving_path):
    os.makedirs(saving_path)
file_error = open('{}/00_error.txt'.format(saving_path),'w')



dRy = []
non_dRy = []

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

    file_variant = open(filepath_variat, 'r')
    variants = file_variant.read().split('\n')
    for variant in variants:
        if not variant.startswith('missense'):
            continue

        v = variant.split('\t')
        generic, before_codon, before_amino, after_codon, after_amino = v[5], v[6], v[7], v[8], v[9]
        ac = int(v[18].split('=')[1])

        if before_amino != 'R':
            continue

        if generic == '3x50':
            dRy.append([[before_codon, after_codon, after_amino],ac])
        else:
            non_dRy.append([[before_codon, after_codon, after_amino],ac])
    
print(dRy)


amino_ls = ['G','W','C','I','M','T','K','S','L','P','Q','H']
amino_colors = ['#ff0101ff',
          '#ff9b2aff',
          '#ffbf07ff',
          '#a7ff58ff',
          '#27ba16ff',
          '#59ffe3ff',
          '#4fc0ffff',
          '#518bffff',
          '#332affff',
          '#7118ffff',
          '#ba04ffff',
          '#ff2abcff',
          ]
codon_ls = ['GGT','GGC','GGA','GGG',
            'TGG',
            'TGT','TGC',
            'ATA',
            'ATG',
            'ACA','ACG',
            'AAA','AAG',
            'AGT','AGC',
            'CTT','CTC','CTA','CTG',
            'CCT','CCC','CCA','CCG',
            'CAA','CAG',
            'CAT','CAC',
            ]

codon_colors = ['#ff0101ff','#ff6060ff','#ffa1a1ff','#ffd7d7ff',
          '#ff9b2aff',
          '#ffbf07ff','#ffe69fff',
          '#a7ff58ff',
          '#27ba16ff',
          '#59ffe3ff','#a9fff1ff',
          '#4fc0ffff','#88d4ffff',
          '#518bffff','#a0c0ffff',
          '#332affff','#6059ffff','#8781ffff','#cfcdffff',
          '#7118ffff','#924effff','#b180ffff','#cdaeffff',
          '#ba04ffff','#dc81ffff',
          '#ff2abcff','#ff7fd7ff'
          ]
print(len(codon_ls))

dRy_count_amino = [0 for _ in range(len(amino_ls))]
non_dRy_count_amino = [0 for _ in range(len(amino_ls))]
dRy_count_codon = [0 for _ in range(len(codon_ls))]
non_dRy_count_codon = [0 for _ in range(len(codon_ls))]

for ls in dRy:
    amino_index = amino_ls.index(ls[0][2])
    dRy_count_amino[amino_index] = dRy_count_amino[amino_index] + ls[1]
    
    codon_index = codon_ls.index(ls[0][1])
    dRy_count_codon[codon_index] = dRy_count_codon[codon_index] + ls[1]

for ls in non_dRy:
    amino_index = amino_ls.index(ls[0][2])
    non_dRy_count_amino[amino_index] = non_dRy_count_amino[amino_index] + ls[1]
    
    codon_index = codon_ls.index(ls[0][1])
    non_dRy_count_codon[codon_index] = non_dRy_count_codon[codon_index] + ls[1]


counts = [dRy_count_amino, non_dRy_count_amino, dRy_count_codon, non_dRy_count_codon]
print(counts)
new_counts = [[],[],[],[]]
for j in range(4):
    sum_ = sum(counts[j])
    for k in range(len(counts[j])):
        new_counts[j].append(round(counts[j][k]/sum_*100, 1))



fig, ax = plt.subplots(1,1, figsize=(8,2))

for n in range(len(counts[0])):
    ax.plot([sum(new_counts[0][:n]), sum(new_counts[1][:n])], [1.6,0.4], '--k', zorder=1, linewidth = 0.5)
ax.plot([100,100],[1.6,0.4], '--k', zorder=1, linewidth = 0.5)

for n in range(len(counts[0])):
    ax.barh(2, new_counts[0][n], left=sum(new_counts[0][:n]), color=amino_colors[n], height = 0.8)
    ax.barh(0, new_counts[1][n], left=sum(new_counts[1][:n]), color=amino_colors[n], height = 0.8)


ax.set_ylim(-1,3)

#tickを消す
ax.set_xticks(np.arange(0))
ax.set_yticks(np.arange(0))
#囲い線を消す
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)


fig.savefig('{}/amino.svg'.format(saving_path))
plt.close()



fig, ax = plt.subplots(1,1, figsize=(8,4))

for n in range(len(counts[2])):
    ax.plot([sum(new_counts[2][:n]), sum(new_counts[3][:n])], [1.6,0.4], '--k', zorder=1, linewidth = 0.5)
ax.plot([100,100],[1.6,0.4], '--k', zorder=1, linewidth = 0.5)

for n in range(len(counts[2])):
    ax.barh(2, new_counts[2][n], left=sum(new_counts[2][:n]), color=codon_colors[n], height = 0.8)
    ax.barh(0, new_counts[3][n], left=sum(new_counts[3][:n]), color=codon_colors[n], height = 0.8)


ax.set_ylim(-1,3)

#tickを消す
ax.set_xticks(np.arange(0))
ax.set_yticks(np.arange(0))
#囲い線を消す
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)


fig.savefig('{}/codon.svg'.format(saving_path))
plt.close()