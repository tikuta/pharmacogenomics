ncbiからgpcrのリストを取得'/data/ncbi/GPCRTargets.csv'
'https://www.guidetopharmacology.org/GRAC/GPCRListForward?class=A'

common
1_0_get_json.py 
    download json file from EMSEMBLE by using gene name list from PubMed
2_0_list_gene_info.py 
    make a list of gene region from the json file
3_0_list_isoforms_info.py 
    make a list of isoforms region from the json file

filter
1_0_filter_vcf_on_GPCR_gene.py
    filter the variants of dataset according to its position
2_0_filter_vcf_on canoncal.py
    filter the variant that locate on the canonical isoforms region
3_0_filter_vcf_on_classA_canonical_isoform.py
    filter the variant that locate on the canonical isoforms region of class A GPCR
    the information of class are gained from GPCR db
    Since the name of GPCRs on the list of GPCR classes is common name, list the link common name and gene name are made by hand

classify
1_0_class_into_gpcr.py
    classify the variants into each GPCRs

list
1_0_list_each_variant.py
    In the vcf files, the variant that locate at same position are written in a same low.
    In this script, the variants are devided into different low so that the number of unique variant can be known.

common
4_0_devide_into_chromosome.py
    Downloaded GRCh38 human genome assembly contains all chromosomes.
    The file is too big and using as it is takes long calculation time.
    To avoid that, this script devides the assembly into each chromosomes.
    There are several genome assembly in this file, but only use 'Primary assembly' starts with 'NC'.

list
2_0_convert_into_amino.py
    
    To use the imformation of generic number, additional calculation has been done in 'data/GPCRdb/generic_number'.

graph
Fig1/1_bar.py

Fig_2/b_bar_sense.py

Fig_2/c_bar_stracture_sense.py
    !!need to check why the total number is different from the past result.


filter
4_0_filter_vcf_on_other_classes_canonical_isoform.py
    To make the graph('Fig_3/a_length_variant), make the vcf files of variants for each class.

classify
2_0_class_into_gpcr_not_classA.py
    classify variants into each GPCR which does not belong to class A.

list
3_0_list_each_variant_not_classA.py
    Copied the script of 1_0_list_each_variant.

4_0_convert_into_amino.py
    Copied the script of 2_0_list_gene_info.

graph
Fig_3/a_length_variant.py
    Scatter figure to show the correlation between protein length and the number of variant.
    Currently, this figure is showing only missense variant.

Fig3/c_top_AC_variants.py
    Bar graph show the missense variant which is dominant in Japanese group(AF>0.5)
    !!GPR42 has problem and not solved yet.

fig_4/a_variant_accumulation.py
    Bar graph shows a proportion of the number of GPCR which has the variation at the specific generic number position.
    The generic numbers are filterd by the threshold (50% of the number of class A GPCRs).
    
fig_4/b-e_dRy_codon.py
    Check how the codons or amino acids of arginine of DRY motif or other positions are changed.

fig_4/f_drug_interaction.py
    This analysis shows the correlation between preference of G protein and variant accumulation at G protein interacting position.
    To determine G protein interacting positions, used this article (https://www.biorxiv.org/content/10.1101/2022.09.24.508774v1).

fig_5/a_compare_1KGP.py
    This scatter graph shows the correlation between 1KGP dataset and 14KJPN dataset.

    