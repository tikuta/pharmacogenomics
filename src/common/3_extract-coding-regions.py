#!/usr/bin/env python3

import os
import json
import glob

saving_path = '../../results/3_canonical-coding-regions'
if not os.path.exists(saving_path):
    os.makedirs(saving_path)

# write on file
file_write = open('{}/canonical-isoforms.txt'.format(saving_path),'w')
for p in glob.glob("../../results/1_genes/*.json"):
    gene = os.path.splitext(os.path.basename(p))[0]

    j = json.load(open(p, 'r'))
    biotype = j['biotype']
    if biotype != "protein_coding": # some pseudogenes should be removed
        print("Gene {} is removed since biotype = {}".format(gene, biotype))
        continue

    chr_num, strand = j['seq_region_name'], j['strand']
    
    transcripts = j['Transcript']
    for t in transcripts:
        is_canonical = t['is_canonical'] # 0 or 1
        if int(is_canonical) == 0: # we consider canonical isoforms only
            continue

        if not 'Translation' in t:
            raise Exception("No translation for " + gene) # this should NOT occur since "canonical" assumes proteins 
        
        enst_id, ensp_id = t['Translation']['Parent'], t['Translation']['id']

        protein_start, protein_end = t['Translation']['start'], t['Translation']['end']

        exon_ranges = []
        exons = t['Exon']
        for e in exons:
            start, end = e['start'], e['end']
            exon_ranges.append('{}:{}'.format(start, end))

        write = [gene, ensp_id, chr_num, str(strand)] + exon_ranges
        write = '\t'.join(write)
        file_write.write(write + '\n')
