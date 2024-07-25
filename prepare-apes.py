#!/usr/bin/env python3
import gpcrdb
from ape import *
import vcf
import os
from config import CHIMPANZEE_BLOCK_LIST

def main():
    for gpcrdb_entry in gpcrdb.get_filtered_receptor_list():
        homology = EnsemblHomology(gpcrdb_entry, GreatApe.chimpanzee)
        for human_gene_id, chimp_gene_id in homology.gene_ids.items():
            print(gpcrdb_entry.entry_name, human_gene_id, "->", chimp_gene_id)
            if chimp_gene_id in CHIMPANZEE_BLOCK_LIST:
                continue
            ensembl_entry = EnsemblHomologGene(chimp_gene_id, GreatApe.chimpanzee, gpcrdb_entry, human_gene_id)

            gagp = os.path.join("data", "GAGP", "crossmap_sorted.vcf.gz")
            cds_vcf = os.path.join(ensembl_entry.dirpath, f"{chimp_gene_id}_CDS.vcf")
            vcf.filter_vcf([gagp], ensembl_entry.ordered_coding_regions, cds_vcf)

            with open(ensembl_entry.annotated_csv_path, 'w') as f:
                f.write(Annotation.header() + '\n')

                for snv in vcf.iterate_vcf(cds_vcf):
                    anno = ensembl_entry.annotate(snv)
                    f.write(anno.to_csv_line() + '\n')

if __name__ == '__main__':
    main()