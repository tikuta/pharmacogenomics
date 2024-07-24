#!/usr/bin/env python3
import gpcrdb
from ape import EnsemblHomology, EnsemblHomologGene, GreatApe
import vcf
import os

def main():
    for gpcrdb_entry in gpcrdb.get_filtered_receptor_list():
        homology = EnsemblHomology(gpcrdb_entry, GreatApe.chimpanzee)
        print(gpcrdb_entry.entry_name)
        for human_gene_id, chimp_gene_id in homology.gene_ids.items():
            ensembl_entry = EnsemblHomologGene(chimp_gene_id, GreatApe.chimpanzee, gpcrdb_entry, human_gene_id)

            gagp = os.path.join("data", "GAGP", "crossmap_sorted.vcf.gz")
            cds_vcf = os.path.join(ensembl_entry.dirpath, f"{chimp_gene_id}_CDS.vcf")
            vcf.filter_vcf([gagp], ensembl_entry.ordered_coding_regions, cds_vcf)

if __name__ == '__main__':
    main()