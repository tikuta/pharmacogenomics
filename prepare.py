#!/usr/bin/env python3

import os
import gpcrdb
import ensembl
import vcf
import alphamissense
from utils import Region
from config import *

def main():
    # You wil need *.vcf.gz.tbi files in the data directory.
    vcfs_jpn = [
        "./data/54KJPN/tommo-54kjpn-20230626-GRCh38-af-autosome.vcf.gz", 
        "./data/54KJPN/tommo-54kjpn-20230626-GRCh38-af-chrX_PAR2.vcf.gz",
        # No GPCR genes in PAR3 vcf
    ]
    vcfs_global = [
        "./data/1KGP/ALL.chr1.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/1KGP/ALL.chr2.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/1KGP/ALL.chr3.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/1KGP/ALL.chr4.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/1KGP/ALL.chr5.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/1KGP/ALL.chr6.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/1KGP/ALL.chr7.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/1KGP/ALL.chr8.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/1KGP/ALL.chr9.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/1KGP/ALL.chr10.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/1KGP/ALL.chr11.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/1KGP/ALL.chr12.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/1KGP/ALL.chr13.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/1KGP/ALL.chr14.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/1KGP/ALL.chr15.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/1KGP/ALL.chr16.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/1KGP/ALL.chr17.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/1KGP/ALL.chr18.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/1KGP/ALL.chr19.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/1KGP/ALL.chr20.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/1KGP/ALL.chr21.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/1KGP/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/1KGP/ALL.chrX.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
    ]

    #alphamissense.extract()

    for gpcrdb_entry in gpcrdb.get_filtered_receptor_list():
        print(gpcrdb_entry)

        ensembl_entry = ensembl.EnsemblGeneEntry(gpcrdb_entry)
        ensembl_entry.visualize()

        gene_region = Region(ensembl_entry.region.chromosome, ensembl_entry.region.start, ensembl_entry.region.end)
        vcf.filter_vcf(vcfs_jpn, [gene_region], os.path.join(gpcrdb_entry.dirpath, VCF_JPN_GENE_FILENAME))
        vcf.filter_vcf(vcfs_global, [gene_region], os.path.join(gpcrdb_entry.dirpath, VCF_GLOBAL_GENE_FILENAME))

        vcf_jpn_path = os.path.join(gpcrdb_entry.dirpath, VCF_JPN_CDS_FILENAME)
        vcf_global_path = os.path.join(gpcrdb_entry.dirpath, VCF_GLOBAL_CDS_FILENAME)
        
        vcf.filter_vcf(vcfs_jpn, ensembl_entry.ordered_coding_regions, vcf_jpn_path)
        vcf.filter_vcf(vcfs_global, ensembl_entry.ordered_coding_regions, vcf_global_path)

        csv_jpn_path = os.path.join(gpcrdb_entry.dirpath, CSV_JPN_CDS_FILENAME)
        csv_global_path = os.path.join(gpcrdb_entry.dirpath, CSV_GLOBAL_CDS_FILENAME)

        with open(csv_jpn_path, 'w') as f:
            f.write(ensembl.Annotation.header() + '\n')

            for snv in vcf.iterate_vcf(vcf_jpn_path, '54KJPN'):
                anno = ensembl_entry.annotate(snv)
                f.write(anno.to_csv_line() + '\n')

        with open(csv_global_path, 'w') as f:
            f.write(ensembl.Annotation.header() + '\n')

            for snv in vcf.iterate_vcf(vcf_global_path, '54KJPN'):
                anno = ensembl_entry.annotate(snv)
                f.write(anno.to_csv_line() + '\n')
"""
            # Annotate variations on CDS
            jpn_csv = open(os.path.join(dpath, CSV_JPN_CDS_FILENAME), 'w')
            global_csv = open(os.path.join(dpath, CSV_GLOBAL_CDS_FILENAME), 'w')
            header = ["Var_Type", "Chr", "Pos", "rsID", "Res_Num"]
            header += ["Ref_Base", "Ref_Codon", "Ref_AA"]
            header += ["Alt_Base", "Alt_Codon", "Alt_AA"]
            header += ["Seg", "Generic_Num"]
            header += ["AC", "AN", "AF", "pathogenicity"]
            header += ["AC_XX", "AN_XX", "AF_XX"]
            header += ["AC_XY", "AN_XY", "AF_XY"]
            jpn_csv.write('#' + '\t'.join(header) + '\n')
            global_csv.write('#' + '\t'.join(header) + '\n')

            jpn_vcf = open(os.path.join(dpath, VCF_JPN_CDS_FILENAME))
            global_vcf = open(os.path.join(dpath, VCF_GLOBAL_CDS_FILENAME))

            for f, csv_file in zip([jpn_vcf, global_vcf], [jpn_csv, global_csv]):
                for l in f.readlines():
                    if len(l) == 0:
                        continue

                    var = Variation(l)
                    if not var.passed:
                        continue
                    for snv in var.snvs():
                        anno = gene.annotate(snv)
                        try:
                            res = matched['residues'][anno.res_num - 1]
                        except:
                            print(anno.res_num, snv.position, snv.rsid, snv.ref, snv.alt)
                        assert(res['ensembl_sequence_number'] == anno.res_num)
                        seg = res['segment']
                        generic_num = res['generic_number']

                        pathogenicity = None
                        if anno.var_type == VariationType.MISSENSE:    
                            for am in ams:
                                if snv.chromosome == am.chromosome and snv.position == am.position and snv.alt == am.alt:
                                    if snv.ref != am.ref:
                                        print("GRCh38 vs AlphaMissense base unmatch!", snv.ref, "vs.", am.ref)
                                    sub = anno.ref_aa + str(anno.res_num) + anno.alt_aa
                                    if sub != am.variant:
                                        print("GRCh38 vs AlphaMissense substitution unmatch!", sub, "vs.", am.variant)
                                    pathogenicity = am.pathogenicity

                            if not pathogenicity:
                                print("No corresponding variant in AlphaMissense!", snv)
                        else:
                            pathogenicity = -1

                        cols = [anno.var_type, snv.chromosome, snv.position, snv.rsid, anno.res_num]
                        cols += [snv.ref, anno.ref_codon, anno.ref_aa]
                        cols += [snv.alt, anno.alt_codon, anno.alt_aa]
                        cols += [seg, generic_num]
                        cols += [snv.AC, snv.AN, snv.AF, pathogenicity]
                        csv_file.write('\t'.join([str(col) for col in cols]) + '\n')
"""
if __name__ == '__main__':
    main()