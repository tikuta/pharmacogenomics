#!/usr/bin/env python3

import os
import json

#declare place to save
saving_path = '../../results/3_0_list_isoforms_info'
if not os.path.exists(saving_path):
    os.makedirs(saving_path)
file_error = open('{}/00_error.txt'.format(saving_path),'w')

#make gpcr genename list
file_gpcr = open('../../data/NCBi/GPCRTargets.csv','r')
lines = file_gpcr.read().split('\n')
gpcrs = []
for line in lines:
    if len(line) <= 1:
        continue
    if line.startswith('"Type"'):
        continue
    l = line.split('","')
    gene_name = l[10]
    if gene_name == '':
        file_error.write('{}:\tno gene name\n'.format(line))
        continue
    gpcrs.append(gene_name)
print(gpcrs)

# write on file
file_write = open('{}/isoforms.txt'.format(saving_path),'w')
for gpcr in gpcrs:
    print(gpcr)
    filepath_dic = '../../results/1_0_get_json_file/{}.json'.format(gpcr)
    if not os.path.exists(filepath_dic):
        file_error.write('{}:\tno json file\n'.format(gpcr))
        continue
    file_dic = open(filepath_dic, 'r')
    DIC = json.load(file_dic)
    chr_num, strand = DIC['seq_region_name'],DIC['strand']
    TRANS = DIC['Transcript']

    for j in range(len(TRANS)):#for each isoform
        if not 'Translation' in TRANS[j]:
            file_error.write('{},Exon{}:\t no translation(protein) range info\n'.format(gpcr, str(j)))
            continue
        enst_id, ensp_id = TRANS[j]['Translation']['Parent'], TRANS[j]['Translation']['id']

        is_canonical = TRANS[j]['is_canonical']
        protein_start, protein_end = TRANS[j]['Translation']['start'], TRANS[j]['Translation']['end']
        exons = []
        EXON = TRANS[j]['Exon']
        for n in range(len(EXON)):#for each exon
            start, end = EXON[n]['start'], EXON[n]['end']
            #check if it includes no CDS region
            if start <= protein_start:
                start = protein_start
                if  end <= protein_start:
                    continue
            if end >= protein_end:
                end = protein_end
                if start >= protein_end:
                    continue
            exon = '{}:{}'.format(str(start), str(end))
            exons.append(exon)

        
        write = [gpcr, ensp_id, chr_num, str(strand), str(j), str(is_canonical)] + exons#"j" indicates the number of isoform in json file
        write = '\t'.join(write)
        file_write.write(write + '\n')
