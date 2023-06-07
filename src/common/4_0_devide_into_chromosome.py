#!/usr/bin/env python3

import os
import gzip

saving_path = '../../results/4_0_devide_into_chromosome'
if not os.path.exists(saving_path):
    os.makedirs(saving_path)

file_genome = gzip.open('../../data/GRCh38/genome_assemblies_genome_fasta/ncbi-genomes-2023-04-27/GCF_000001405.26_GRCh38_genomic.fna.gz','rt')
line = file_genome.readline()

file_index = open('{}/index.txt'.format(saving_path),'w')

n = 1
tf = False
while line:
    if line.startswith('>'):
        file_index.write(line)
        tf = False
    if line.startswith('>NC'):
        chr_num = str(n)
        n = n + 1
        if chr_num == '23':
            chr_num = 'X'
        elif chr_num == '24':
            chr_num = 'Y'
        elif chr_num == '25':
            chr_num = '_mitochondrion'
        file_write = open('{}/chr{}.txt'.format(saving_path, chr_num), 'w')
        print(line)
        tf = True
    if tf == True:
        file_write.write(line)
    line = file_genome.readline()
    continue