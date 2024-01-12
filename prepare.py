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
        vcf.filter_vcf(vcfs_jpn, [gene_region], gpcrdb_entry.japan_gene_vcf_path, force=True)
        vcf.filter_vcf(vcfs_global, [gene_region], gpcrdb_entry.global_gene_vcf_path, force=True)
        
        vcf.filter_vcf(vcfs_jpn, ensembl_entry.ordered_coding_regions, gpcrdb_entry.japan_cds_vcf_path, force=True)
        vcf.filter_vcf(vcfs_global, ensembl_entry.ordered_coding_regions, gpcrdb_entry.global_cds_vcf_path, force=True)

        with open(gpcrdb_entry.japan_cds_csv_path, 'w') as f:
            f.write(ensembl.Annotation.header() + '\n')

            for snv in vcf.iterate_vcf(gpcrdb_entry.japan_cds_vcf_path, '54KJPN'):
                anno = ensembl_entry.annotate(snv)
                f.write(anno.to_csv_line() + '\n')

        with open(gpcrdb_entry.global_cds_csv_path, 'w') as f:
            f.write(ensembl.Annotation.header() + '\n')

            for snv in vcf.iterate_vcf(gpcrdb_entry.global_cds_vcf_path, '1KGP'):
                anno = ensembl_entry.annotate(snv)
                f.write(anno.to_csv_line() + '\n')

if __name__ == '__main__':
    main()