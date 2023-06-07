#!/usr/bin/env python3

#missense のみになってない

import os
import numpy as np
import matplotlib.pyplot as plt

saving_path = '../../../results/graph/Fig_4/f_drug_interaction'
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

print(GENERIC)

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


#---------------------------------------------------------------------------------------------------------
# copied from '/Volumes/T7/final/genome_variation/38kjpn/src/graph/38kjpn_special/O1_Gprotein_coupling_classify.py'
file = open('../../../data/GPCRdb/G_coupling/families_coupling.csv', 'r')
lines = file.read().split('\n')

file_g = open('../../../data/GPCRdb/G_coupling/handmade_gene_name_Gcoupling.txt','r')
names = file_g.read().split('\n')
NAME_G = {}
for name in names:
    n = name.split(',')
    NAME_G[n[0].replace('\ufeff','')] = NAME_G.get(n[0].replace('\ufeff',''),'')
    NAME_G[n[0].replace('\ufeff','')] = n[1] 

DICS = {'Bouvier':{}, 'GPCRdb':{}, 'Inoue':{}}
for line in lines:
    if len(line) == 0:
        continue
    l = line.split(',')
    data_name, c, gene_name, nums = l[1], l[3], l[5], l[9:13]
    if gene_name in NAME_G:
        gene_name = NAME_G[gene_name]
    if not c == 'A':
        continue
    if not data_name in DICS:
        continue
    DICS[data_name][gene_name] = nums
#print(DICS)

#do not untag this section(hand made list will be deleted)
#file = open('../../../data/GPCRdb/G_coupling/handmade_gene_name_Gcoupling.txt','w')
#for gene_ in DICS['Inoue']:
    #file.write(gene_ + ',\n')

#print(DICS['Inoue'])


#------------------------------------------------------------------------------------------------------

target_generics = ['2x37','2x39','3x50','3x53','3x54','34x50','34x51','34x54','34x55','5x65','5x68','6x29','6x32','6x33','6x36','6x37','7x56','8x47','8x48','8x49']

COUNT = {}
for target in target_generics:
    COUNT[target] = COUNT.get(target, [[0,0,0,0],[0,0,0,0]])
count_gpcrs = 0

check = []

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

    filepath_variant = '../../../results/list/2_0_convert_into_amino/38KJPN/{}.txt'.format(gene_name)
    if not os.path.exists(filepath_variant):
        file_error.write('{}:\t no variant file exists (may be not in class A)\n'.format(gene_name))
        continue


    #change gene_name to suitable name for the DICS
    if not gene_name in DICS['Inoue']:
        file_error.write('{}:\tnot g_coupling list\n'.format(gene_name))
        if gene_name.startswith('GPR'):
            gene_name = gene_name.replace('R','')
            if not gene_name in DICS['Inoue']:
                continue
        else:
            continue
    check.append(gene_name)# Use this check list to check if any gpcr that is not counted even though it exists on the dictionary.
    
    coupling = DICS['Inoue'][gene_name]

    orginal_gene_name = i[0]

    # check if the gpcr has the generic number at specific position. If it has, "COUNT" will be added.
    GENERIC_EACH = GENERIC[orginal_gene_name]
    generics = np.array(list(GENERIC_EACH.values()))
    generics_tp = np.transpose(generics)[0]
    for target in target_generics:
        if not target in generics_tp:
            continue
        for n in range(4):
            if coupling[n] == "1'":
                COUNT[target][0][n] = COUNT[target][0][n] + 1


    file_variant = open(filepath_variant, 'r')
    variants = file_variant.read().split('\n')

    for variant in variants:
        if len(variant) == 0:
            continue
        v = variant.split('\t')
        generic = v[5]
        if not v[5] in target_generics:
            continue
        for n in range(4):
            if coupling[n] == "1'":
                COUNT[generic][1][n] = COUNT[generic][1][n] + 1

    
print(COUNT)


# To check if the all GPCR in the dictionary are used.------------------------------------------------    
for gene_name in DICS['Inoue']:
    if gene_name not in check:
        file_error.write('{}:\t not used \n'.format(gene_name))
    

#---------------------------------------------------------------------------------------------------------
keys = list(COUNT.keys())
values_ = np.array(list(COUNT.values()))
print(values)
values = []
for value_ in values_:
    value = []
    for n in range(4):
        value.append(round(value_[1][n] / value_[0][n] * 100, 1))
    values.append(value)
print(values)

labels = ['Gs','Gi','Gq','G12']

angles = np.linspace(0, 2*np.pi, len(labels)+1, endpoint=True)
rgrids = [0,20,40,60,80,100]

fig = plt.figure(facecolor='w')

ax = fig.add_subplot(1, 1, 1, polar = True)



# 多角形の目盛線を引く
for grid_value in rgrids:
    grid_values = [grid_value] * (len(labels)+1)
    ax.plot(angles, grid_values, color="gray",  linewidth=0.5)

for t in rgrids:
    ax.text(x=0, y=t, s=t)
# rの範囲を指定

for value in values:
    radar_values = np.concatenate([value, [value[0]]])
    ax.plot(angles, radar_values)



ax.set_rlim([min(rgrids), max(rgrids)])


ax.set_thetagrids(angles[:-1]*180/np.pi, labels)
ax.set_rgrids([])
ax.spines['polar'].set_visible(False)
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)

plt.close()

# make the image for each position-------------------------------------------------------------------


for value in values:
    fig = plt.figure(facecolor='w')

    ax = fig.add_subplot(1, 1, 1, polar = True)



    # 多角形の目盛線を引く
    for grid_value in rgrids:
        grid_values = [grid_value] * (len(labels)+1)
        ax.plot(angles, grid_values, color="gray",  linewidth=0.5)

    for t in rgrids:
        ax.text(x=0, y=t, s=t)
    # rの範囲を指定


    radar_values = np.concatenate([value, [value[0]]])
    ax.plot(angles, radar_values)



    ax.set_rlim([min(rgrids), max(rgrids)])


    ax.set_thetagrids(angles[:-1]*180/np.pi, labels)
    ax.set_rgrids([])
    ax.spines['polar'].set_visible(False)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_title("{}".format(keys[values.index(value)]), pad=20)
    fig.savefig('{}/each_{}.svg'.format(saving_path, keys[values.index(value)]))

    plt.close()
