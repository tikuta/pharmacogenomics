#!usr/bin/env python3

import os
import numpy as np
import sys
import matplotlib.pyplot as plt
import gzip

saving_path = '../../../results/graph/Fig_1/1_bar'
if not os.path.exists(saving_path):
    os.makedirs(saving_path)
file_error = open('{}/00_error.txt'.format(saving_path),'w')

filepath_all = ['../../../data/38KJPN/tommo-38kjpn-20220929-GRCh38-gf-chrX_PAR2.vcf.gz','../../../data/38KJPN/tommo-38kjpn-20220929-GRCh38-gf-chrX_PAR3.vcf.gz','../../../data/38KJPN/tommo-38kjpn-20220630-GRCh38-gf-autosome.vcf.gz']
filepath_gene = '../../../results/filter/1_0_fileter_vcf_on_GPCR_gene/38KJPN/on_GPCR_gene.vcf'
filepath_cds = '../../../results/filter/2_0_fileter_vcf_on_canonical_isoform/38KJPN/on_canonical_isoforms.vcf'




# total, SNVs, indel
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

tf = True
if tf == True:
    file = open(filepath_gene,'r')
    line = file.readline()
    gene_c = [0,0,0]
    while line:
        if line.startswith('#'):
            line = file.readline()
            continue
        count(line, gene_c)
        line = file.readline()
    print(gene_c)

    file = open(filepath_cds,'r')
    line = file.readline()
    cds_c = [0,0,0]
    while line:
        if line.startswith('#'):
            line = file.readline()
            continue
        count(line, cds_c)
        line = file.readline()
    print(cds_c)


    all_c = [0,0,0]
    for filepath in filepath_all:
        file = gzip.open(filepath, 'rt')
        line = file.readline()
        while line:
            if line.startswith('#'):
                line = file.readline()
                continue
            count(line, all_c)
            line = file.readline()
    print(all_c)


file_row = open('{}/row_numbers.csv'.format(saving_path),'w')
file_row.write('all:\t{}\n'.format('\t'.join(all_c)))
file_row.write('gene:\t{}\n'.format('\t'.join(gene_c)))
file_row.write('CDS:\t{}\n'.format('\t'.join(cds_c)))

file_statistic = open('{}/statistic.csv'.format(saving_path), 'w')
file_statistic.write('.\tother than GPCR genome\tGPCR genome\tother than CDS region\tCDS region\n')


#all_c = [207435063,207435063,207435063]#仮置き
#gene_c =[1437484, 1346522, 222552]
#cds_c = [37744, 37097, 1791]

names = ['total', 'SNVs', 'Indel']
for n in range(3):
    file_statistic.write(names[n] + '\t')
    fig, ax = plt.subplots(1,1,figsize=(4,2))


    xs_0 = [cds_c[n]/gene_c[n]*100,(gene_c[n]-cds_c[n])/gene_c[n]*100]
    xs_1 = [gene_c[n]/all_c[n]*100, (all_c[n] - gene_c[n])/all_c[n]*100]

    #区分線を引く
    ind = np.arange(2)
    width = 0.8
    ind_p = ind + width/2
    ind_m = ind - width/2
    ind_line = np.sort(np.concatenate([ind_p, ind_m]))[::-1]
    A_line = np.insert(xs_0, np.arange(2), xs_0)
    ax.plot(A_line, ind_line, '--k', zorder=1)
    ax.plot([0,0], [1-width, 0+width], '--k', zorder=1)
    

    #CDSに存在するバリアントのサイト数を示す
    colors = ['tomato','lightgrey']
    for i in range(2):
        p = ax.barh(0, round(xs_0[i],1), left=sum(xs_0[:i]), color = colors[i], height = width)
        file_statistic.write(str(xs_0[i]) + '\t')

    labels = [ 'CDS region\n{}%'.format(round(xs_0[0],1)),'\n{}%'.format(round(xs_0[1],1))]
    ax.text(xs_0[0]/2, 0, labels[0], ha='center', va='center', fontfamily='Arial')
    ax.text(xs_0[1]/2+xs_0[0], 0, labels[1], ha='center', va='center', fontfamily='Arial')


    #GPCRの遺伝子に存在するバリアントのサイト数を示す
    colors = ['red','lightgrey']
    for i in range(2):
        p = ax.barh(1, round(xs_1[i],1), left=sum(xs_1[:i]), color = colors[i], height=width)
        file_statistic.write(str(xs_1[i]) + '\t')

    labels = ['GPCR genome\n{}%'.format(round(xs_1[0],1)),'\n{}%'.format(round(xs_1[1],1))]
    ax.text(xs_1[0]/2, 1, labels[0], ha='center', va='center', fontfamily='Arial')
    ax.text(xs_1[1]/2+xs_1[0], 1, labels[1], ha='center', va='center', fontfamily='Arial')



    #tickを消す
    ax.set_xticks(np.arange(0))
    ax.set_yticks(np.arange(0))
    #囲い線を消す
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    fig.savefig('{}/cds_barh_{}.svg'.format(saving_path, names[n]), dpi=300)

    file_statistic.write('\n')
