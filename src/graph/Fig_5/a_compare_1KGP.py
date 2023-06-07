#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt

saving_path = '../../../results/graph/Fig_5/a_compare_1KGP'
if not os.path.exists(saving_path):
    os.makedirs(saving_path)
file_error = open('{}/00_error.txt'.format(saving_path),'w')

COUNT = {}

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

    filepath_variant_38kjpn = '../../../results/list/2_0_convert_into_amino/38KJPN/{}.txt'.format(gene_name)
    filepath_variant_1kgp = '../../../results/list/2_0_convert_into_amino/1KGP/{}.txt'.format(gene_name)

    if not os.path.exists(filepath_variant_38kjpn):
        file_error.write('{}:\tno 38KJPN variant data\n'.format(gene_name))
        continue
    if not os.path.exists(filepath_variant_1kgp):
        file_error.write('{}:\tno 1KGP variant data\n'.format(gene_name))
        continue

    file_variant_38kjpn = open(filepath_variant_38kjpn,'r')
    variants_38kjpn = file_variant_38kjpn.read().split('\n')
    file_variant_1kgp = open(filepath_variant_1kgp, 'r')
    variants_1kgp = file_variant_1kgp.read().split('\n')

    COUNT_EACH = {}
    
    for variant in variants_38kjpn:
        if len(variant) == 0:
            continue
        if not variant.startswith('missense'):
            continue
        v = variant.split('\t')
        key = '{}_{}{}{}_{}'.format(gene_name, v[7], v[2], v[9], v[5])
        print(key)
        for info in v[10:]:
            if info.startswith('AC'):
                ac = int(info.split('=')[1])
            if info.startswith('AN'):
                an = int(info.split('=')[1])
        value = ac / an
        COUNT_EACH[key] = COUNT_EACH.get(key, [0,0])
        COUNT_EACH[key][0] = value

    for variant in variants_1kgp:
        if not variant.startswith('missense'):
            continue
        v = variant.split('\t')
        key = '{}_{}{}{}_{}'.format(gene_name, v[7], v[2], v[9], v[5])
        if not key in COUNT_EACH:
            continue
        for info in v[18].split(';'):
            if info.startswith('AC'):
                ac = int(info.split('=')[1])
            if info.startswith('AN'):
                an = int(info.split('=')[1])
        value = ac / an
        COUNT_EACH[key][1] = value

    COUNT.update(COUNT_EACH)
    
print(COUNT)
exists_ls = []
nonexists_ls = []

for key in COUNT:
    if COUNT[key][1] == 0:
        COUNT[key].append('inf')
        nonexists_ls.append(COUNT[key])
    else:
        COUNT[key].append(COUNT[key][0]/COUNT[key][1])
        exists_ls.append(COUNT[key])
print(exists_ls)


exists_tp = np.transpose(np.array(exists_ls))
nonexists_tp = np.transpose(np.array(nonexists_ls))


fig = plt.figure(linewidth=1)
grid = plt.GridSpec(6,1)
ax1 = fig.add_subplot(grid[0,0])
ax2 = fig.add_subplot(grid[1:,0])

ax1.scatter(nonexists_tp[0],[0 for _ in range(len(nonexists_tp[0]))], s=1)
ax2.scatter(exists_tp[0], exists_tp[2], s = 1)

ax1.set_yticks(np.arange(1))
ax1.set_yticklabels(['inf.'])
ax1.set_xticks(np.arange(1))

limit = 1
#ax.set_xlim(0,limit)
#ax.set_ylim(0,limit)
ax2.set_xlabel('AF of 38KJPN')
ax2.set_ylabel('AF(38KJPN) / AF(1KGP)')


fig.savefig('{}/relative_value.svg'.format(saving_path))