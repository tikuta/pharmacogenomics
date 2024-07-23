#!/usr/bin/env python3
import gpcrdb
from ape import EnsemblHomology, EnsemblChimpanzeeGene, GreatApe
import vcf
import os

def main():
    for gpcrdb_entry in gpcrdb.get_filtered_receptor_list():
        print(gpcrdb_entry.entry_name)
        homology = EnsemblHomology(gpcrdb_entry)
        for gene_id in homology.gene_ids[GreatApe.chimpanzee]:
            ensembl_entry = EnsemblChimpanzeeGene(gene_id)

            gagp = os.path.join("data", "GAGP", "crossmap_sorted.vcf.gz")
            cds_vcf = os.path.join(ensembl_entry.dirpath, f"{gene_id}_CDS.vcf")
            vcf.filter_vcf([gagp], ensembl_entry.ordered_coding_regions, cds_vcf)


if __name__ == '__main__':
    main()