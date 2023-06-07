#!/usr/bin/env python3
import os

saving_path ='../../results/classify/1_0_class_into_gpcr'
if not os.path.exists(saving_path):
    os.makedirs(saving_path)
file_error = open('{}/00_error.txt'.format(saving_path),'w')

#make dictionary (key = common name, value = gene name)
file_name = open('../../data/GPCRdb/handmade_gene_name.csv','r')
names = file_name.read().split('\n')
NAME = {}
for name in names:
    n = name.split(',')
    NAME[n[0].replace('\ufeff','')] = NAME.get(n[0].replace('\ufeff',''),'')
    NAME[n[0].replace('\ufeff','')] = n[1] 
#print(NAME)

#list up the gpcrs that belong to class A
file_class = open('../../data/GPCRdb/GPCRdb_alignment_A.csv','r')
lines = file_class.read().split('\n')
a_list = []
A_NAME = {}
for line in lines:
    if len(line) == 0:
        continue
    if not line.startswith('[Human]'):
        continue
    receptor = line.split(',')[0].replace('[Human]','').replace('receptor','').replace(' ','')
    if receptor in NAME:
        receptor = NAME[receptor]
    A_NAME[receptor] = A_NAME.get(receptor, False)
    a_list.append(receptor)
#print(a_list)
#print(len(a_list))

# make a dictionary (key = {chrmosome number}_{position on the chromosome}, value = [list of the gpcr gene which locates on the position]})
# used a list for value in case that two gene (plus or minus read) exist at same position of the chromosome
file_gene = open('../../results/3_0_list_isoforms_info/isoforms.txt','r')
genes = file_gene.read().split('\n')
DIC = {}
for gene in genes:
    if len(gene) == 0 :
        continue
    is_canonical = gene.split('\t')[5]
    if is_canonical != '1':
        continue
    gene_name, emsemble_id, chr_m = gene.split('\t')[0],gene.split('\t')[1],gene.split('\t')[2]
    if not gene_name in a_list:
        continue
    A_NAME[gene_name] = True
    filepath = '{}/{}_{}.txt'.format(saving_path, gene_name, emsemble_id)
    regions = gene.split('\t')[6:]
    #print(gs)
    region_list = []
    for region in regions:
        rs = region.split(':')
        ls = []
        for n in range(int(rs[0]),int(rs[1])):
            key = '{}_{}'.format(chr_m, str(n))
            DIC[key] = DIC.get(key, [])
            DIC[key].append(gene_name)
#print(DIC)

datasets = ['1KGP','38KJPN']
# print error gpcr which has not considered (mostly pseudogene?)
# make a filepath list to save
FILE = {}
for common_name in A_NAME:
    if A_NAME[common_name] == False:
        file_error.write(common_name + ':\thas not been considered, because the gene name is not in "../../results/3_0_list_isoforms_info.py/isoforms.txt" \n')
        continue
    FILE[common_name] = FILE.get(common_name, {})
    for dataset in datasets:
        FILE_DATA = FILE[common_name]
        file_path_write = '{}/{}/{}.txt'.format(saving_path,dataset,common_name)
        os.remove(file_path_write) #To avoid writing same variant several times in case such as runnning this script several times, delete the saving file before writing with paramater 'a'.
        file_write = open(file_path_write, 'w')
        FILE_DATA[dataset] = FILE_DATA.get(dataset, '')
        FILE_DATA[dataset] = file_path_write
print(FILE)

# classify the variants
# check if the variant position exists in DIC(GPCR gene position list)
for dataset in datasets:
    filepath = '{}/{}'.format(saving_path,dataset)
    if not os.path.exists(filepath):
        os.makedirs(filepath)
    file_variants = open('../../results/filter/3_0_filter_vcf_on_classA_canonical_isoform/{}/on_classA_canonical_isoforms.vcf'.format(dataset),'r')
    variants = file_variants.read().split('\n')
    for variant in variants:
        if len(variant) == 0:
            continue
        v = variant.split('\t')
        v_chr_m, v_pos = v[0].replace('chr',''), v[1]
        key = '{}_{}'.format(v_chr_m, v_pos)
        if not key in DIC:
            continue
        values = DIC[key]
        for value in values:
            file_path_write = FILE[value][dataset]
            print(file_path_write)
            file = open(file_path_write,'a')#used 'a' parameter for writing in this script.
            file.write(variant +'\n')


# if you open the all saving first at beginning, it will respond error because too many files are opened. This is the reason why the paramater 'a' is used.
