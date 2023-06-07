#!/usr/bin/env python3


import os
import numpy as np
import sys

saving_path = '../../results/list/4_0_convert_into_amino_not_classA'
if not os.path.exists(saving_path):
    os.makedirs(saving_path)
file_error = open('{}/00_error.txt'.format(saving_path),'w')


#----DICTIONARIES--------------------------------------------------------------------------------------------------------

FOR_MINUS = {
    'A':'T', 'T':'A', 'G':'C', 'C':'G',
    'a':'T', 't':'A', 'g':'C', 'c':'G'
}
CAPITALIZE = {
    'a':'A', 't':'T', 'g':'G', 'c':'C',
    'A':'A', 'T':'T', 'G':'G', 'C':'C'
}

DNA = {
        'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C',
        'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C',
        'TTA' : 'L', 'TCA' : 'S', 'TAA' : '*', 'TGA' : '*',
        'TTG' : 'L', 'TCG' : 'S', 'TAG' : '*', 'TGG' : 'W',

        'CTT' : 'L', 'CCT' : 'P', 'CAT' : 'H', 'CGT' : 'R',
        'CTC' : 'L', 'CCC' : 'P', 'CAC' : 'H', 'CGC' : 'R',
        'CTA' : 'L', 'CCA' : 'P', 'CAA' : 'Q', 'CGA' : 'R',
        'CTG' : 'L', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R',

        'ATT' : 'I', 'ACT' : 'T', 'AAT' : 'N', 'AGT' : 'S',
        'ATC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S',
        'ATA' : 'I', 'ACA' : 'T', 'AAA' : 'K', 'AGA' : 'R',
        'ATG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'R',

        'GTT' : 'V', 'GCT' : 'A', 'GAT' : 'D', 'GGT' : 'G',
        'GTC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G',
        'GTA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G',
        'GTG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G'
}

#make dictionary (key = common name, value = gene name)
file_name = open('../../data/GPCRdb/handmade_gene_name_2.csv','r')
names = file_name.read().split('\n')
NAME = {}
for name in names:
    n = name.split(',')
    NAME[n[0].replace('\ufeff','')] = NAME.get(n[0].replace('\ufeff',''),'')
    NAME[n[0].replace('\ufeff','')] = n[1] 

'''
GENERIC = {}
filepath_generics = ['../../data/GPCRdb/generic_number/residue_table_startswithO.csv','../../data/GPCRdb/generic_number/residue_table_startswithP.csv','../../data/GPCRdb/generic_number/residue_table_startswithQ.csv']
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
keys = GENERIC.keys()
print(keys)
print(len(keys))
'''
#----FUNCTIONS------------------------------------------------------------------------------------------------------------



def amino_frame(num):
    '''
    Return the position of a residue, the positon of an amino acid, and frame number
    If the codon were ATG, A = frame 0, T = frame 1, G = frame 2.
    '''
    
    if num % 3 == 1:
        frame = 0
        pos = num // 3 + 1
    if num % 3 == 2:
        frame = 1
        pos = num // 3 + 1
    if num % 3 == 0:
        frame = 2
        pos = num // 3
    return [num, pos, frame]


def get_codon(v_pos):
    '''
    By using the dictionary contains the position and residue info, this function allows to get the codon.
    Using the list 'residue_nubmers', this functiona can return appropriate codon even if the reading direction is minus.
    '''
    info = DIC[v_pos]
    frame = int(info[2])
    r = residue_numbers.index(v_pos)
    if info[2] == 0:
        first = DIC[residue_numbers[r]][3]
        second = DIC[residue_numbers[r+1]][3]
        third = DIC[residue_numbers[r+2]][3]
    if info[2] == 1:
        first = DIC[residue_numbers[r-1]][3]
        second = DIC[residue_numbers[r]][3]
        third = DIC[residue_numbers[r+1]][3]
    if info[2] == 2:
        first = DIC[residue_numbers[r-2]][3]
        second = DIC[residue_numbers[r-1]][3]
        third = DIC[residue_numbers[r]][3]
    #print(first, second, third)
    codon = CAPITALIZE[first] + CAPITALIZE[second] + CAPITALIZE[third]
    return codon, frame


#------------------------------------------------------------------------------------------------------------------------------

file_isoforms = open('../../results/3_0_list_isoforms_info/isoforms.txt', 'r')
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

    #------------------------------------------------------------------------------------------
    # Make a dictionary of conserved receptor
    DIC = {}#key = position in genome, value = [the position of a residue, the positon of an amino acid, and frame number, reisude]
    if m_p == '1':
        n = 1
        for position in positions:
            pos = position.split(':')
            for k in range(int(pos[0]), int(pos[1])+1):
                key = str(k)
                DIC[key] = DIC.get(key, [])
                value = amino_frame(n)
                DIC[key] = value
                n = n + 1

    elif m_p == '-1':
        n = 1
        #positions.reverse()
        for position in positions:
            pos = position.split(':')
            p_nums = list(range(int(pos[0]),int(pos[1])+1))
            p_nums.reverse()
            for p in p_nums:
                key = str(p)
                DIC[key] = DIC.get(key, [])
                value = amino_frame(n)
                DIC[key] = value
                n = n + 1
    #print(DIC)
    residue_numbers = list(DIC.keys())
    int_residue_numbers = []
    for res in residue_numbers:
        int_residue_numbers.append(int(res))
    minimum_pos = int(min(residue_numbers))
    maximum_pos = int(max(residue_numbers))

    #print(residue_numbers)

    #------------------------------------------------------------------------------------------
    # Get a referrence residue of the GPCR CDS.
    # file_genome shows GRCh38 primary genome assembly. It has 80 residues in a low.
    file_genome = open('../../results/4_0_devide_into_chromosome/chr{}.txt'.format(chr_num),'r')
    line = file_genome.readline()

    while line:
        if line.startswith('>'):
            line = file_genome.readline()
            n = 1
            continue
        if n * 80 < minimum_pos:# その行の最大のポジション番号が遺伝子の最初より小さい
            line = file_genome.readline()
            n = n + 1
            continue
        if (n - 1) * 80 > maximum_pos:# その行の最小のポジション番号が遺伝子の最後より大きい
            line = file_genome.readline()
            n = n + 1
            continue

        l = list(line.replace('\n',''))
        for r in range(len(l)):
            residue_number = 80 * (n-1) + r + 1
            if not str(residue_number) in residue_numbers:
                continue
            if m_p == '1':
                DIC[str(residue_number)].append(l[r])
            elif m_p == '-1':
                DIC[str(residue_number)].append(FOR_MINUS[l[r]])
        line = file_genome.readline()
        n = n + 1
    
                    
    #------------------------------------------------------------------------------------------
    # See if the variant locate on the GPCR CDS reagion.
    # Also, verify the codon and the amino acids.
    datasets = ['1KGP','38KJPN']
    for dataset in datasets:
        if not os.path.exists('{}/{}'.format(saving_path, dataset)):
            os.makedirs('{}/{}'.format(saving_path, dataset))

        filepath_variant = '../../results/list/3_0_list_each_variant_not_classA/{}/{}.txt'.format(dataset, gene_name)
        if not os.path.exists(filepath_variant):
            continue
        file_variant = open(filepath_variant, 'r')
        variants = file_variant.read().split('\n')

        file_write = open('{}/{}/{}.txt'.format(saving_path, dataset, gene_name),'w')

        for variant in variants:
            if len(variant) == 0:
                continue
            v = variant.split('\t')


            v_pos = v[1]
            if not v_pos in DIC:
                file_error.write('{},{}:\t not in canonical reagion\n'.format(gene_name, v_pos))
                continue

            if len(v[3]) > len(v[4]):
                file_write.write('deletion\t')
                info = ['.' for _ in range(7)]
                info = '\t'.join(info)
                file_write.write(info + '\t||\t' + variant + '\n')
                continue
            
            if len(v[3]) < len(v[4]):
                file_write.write('insertion\t')
                info = ['.' for _ in range(7)]
                info = '\t'.join(info)
                file_write.write(info + '\t||\t' + variant + '\n')
                continue
            
            if len(v[3]) > 1:
                before = list(v[3])
                after = list(v[4])
                file_error.write('{}:\t longer than 1 ({})\n'.format(gene_name, dataset))
                continue


            codon, frame = get_codon(v_pos)
            before_res, after_res = v[3], v[4]
            if m_p == '-1':
                before_res = FOR_MINUS[before_res]
                after_res = FOR_MINUS[after_res]
            if codon[frame] != before_res:
                file_error.write("{},{}:\t codon doesn't match\n".format(gene_name, v_pos)) 
                continue  
            before_codon = codon 
            before_amino = DNA[before_codon]

            codon_ls = list(codon)
            codon_ls[frame] = after_res
            after_codon = ''.join(codon_ls)
            after_amino = DNA[after_codon]
            '''
            if not gene_name in GENERIC:
                file_error.write('{}:\t not in generic number list \n'.format(gene_name))
                continue
            if DIC[v_pos][1] in GENERIC[gene_name]:
                generic = GENERIC[gene_name][DIC[v_pos][1]][0]
            else:
                generic = '.'
            '''
            generic = '.'

            #convert all info into 'str'
            info = []
            for v_po in DIC[v_pos]:
                if type(v_po) != str:
                    info.append(str(v_po))
                else:
                    info.append(v_po)

            if after_amino == '*':
                file_write.write('nonsense\t')
            elif before_amino == after_amino:
                file_write.write('silent\t')
            else:
                file_write.write('missense\t')

            write = info + [generic, before_codon, before_amino, after_codon, after_amino]
            write = '\t'.join(write)
            file_write.write(write + '\t||\t' + variant+'\n')
            