#!/usr/bin/env python3

import os
import gpcrdb
import uniprot
import ensembl
import vcf
from utils import Region, Gene, Variation
from misc import *

def main():
    vcfs_jpn = [
        "./data/tommo-54kjpn-20230626-GRCh38-af-autosome.vcf.gz", 
        "./data/tommo-54kjpn-20230626-GRCh38-af-chrX_PAR2.vcf.gz",
        # No GPCR gene regions in PAR3 vcf
    ]
    vcfs_global = [
        "./data/ALL.chr1.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/ALL.chr2.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/ALL.chr3.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/ALL.chr4.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/ALL.chr5.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/ALL.chr6.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/ALL.chr7.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/ALL.chr8.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/ALL.chr9.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/ALL.chr10.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/ALL.chr11.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/ALL.chr12.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/ALL.chr13.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/ALL.chr14.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/ALL.chr15.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/ALL.chr16.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/ALL.chr17.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/ALL.chr18.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/ALL.chr19.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/ALL.chr20.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/ALL.chr21.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
        "./data/ALL.chrX.shapeit2_integrated_v1a.GRCh38.20181129.GRCh38.phased.vcf.gz",
    ]

    for d in gpcrdb.get_filtered_receptor_list("receptors.json"):
        entry_name = d['entry_name']
        receptor_class = d['receptor_class']
        dpath = os.path.join(receptor_class, entry_name)
        os.makedirs(dpath, exist_ok=True)

        print(entry_name)

        generic_numbers = gpcrdb.get_generic_number(entry_name, os.path.join(dpath, "gpcrdb.json"))

        accession = d['accession']

        entry = uniprot.get_entry(accession, os.path.join(dpath, "uniprot.json"))
        gene_id = uniprot.uniprot2ensembl(entry)
        gene_info = ensembl.get_gene(gene_id, os.path.join(dpath, "ensembl.json"))

        chromosome = gene_info['seq_region_name']
        gene_start = gene_info['start']
        gene_end = gene_info['end']
        strand = gene_info['strand']
        assert(gene_info['assembly_name'] == 'GRCh38')
        gene_region = Region(chromosome, gene_start, gene_end)
        vcf.filter_vcf(vcfs_jpn, [gene_region], os.path.join(dpath, VCF_JPN_GENE_FILENAME))
        vcf.filter_vcf(vcfs_global, [gene_region], os.path.join(dpath, VCF_GLOBAL_GENE_FILENAME))

        for transcript in gene_info['Transcript']:
            if transcript['is_canonical'] != 1:
                continue

            seqs = ensembl.get_sequence(gene_region, strand, os.path.join(dpath, "sequence.json"))
            coding_regions = ensembl.get_coding_regions(transcript)
            coding_regions.sort(key=lambda r: r.start)

            # Remove STOP codon
            if strand == 1:
                assert(len(coding_regions[-1]) > 3)
                last_region = coding_regions.pop(-1)
                coding_regions.append(Region(last_region.chromosome, last_region.start, last_region.end - 3))
            else:
                assert(len(coding_regions[0]) > 3)
                last_region = coding_regions.pop(0)
                coding_regions.insert(0, Region(last_region.chromosome, last_region.start + 3, last_region.end))

            gene = Gene(gene_region, seqs[0]['seq'], strand, coding_regions)
            gene.visualize(os.path.join(dpath, "gene.png"))
            gene.save_cds(os.path.join(dpath, "cds.json"))

            vcf.filter_vcf(vcfs_jpn, coding_regions, os.path.join(dpath, VCF_JPN_CDS_FILENAME))
            vcf.filter_vcf(vcfs_global, coding_regions, os.path.join(dpath, VCF_GLOBAL_CDS_FILENAME))

            ref_seq = gene.translate()
            assert(ref_seq[0] == 'M')
            matched = gpcrdb.match(ref_seq, generic_numbers, receptor_class, os.path.join(dpath, "match.json"), force=True, alignment_for_human=os.path.join(dpath, "match.txt"))

            # Annotate variations on CDS
            jpn_csv = open(os.path.join(dpath, CSV_JPN_CDS_FILENAME), 'w')
            global_csv = open(os.path.join(dpath, CSV_GLOBAL_CDS_FILENAME), 'w')
            header = ["Var_Type", "Chr", "Pos", "rsID", "Res_Num"]
            header += ["Ref_Base", "Ref_Codon", "Ref_AA"]
            header += ["Alt_Base", "Alt_Codon", "Alt_AA"]
            header += ["Seg", "Generic_Num"]
            header += ["AC", "AN", "AF"]
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
                        cols = [anno.var_type, snv.chromosome, snv.position, snv.rsid, anno.res_num]
                        cols += [snv.ref, anno.ref_codon, anno.ref_aa]
                        cols += [snv.alt, anno.alt_codon, anno.alt_aa]
                        cols += [seg, generic_num]
                        cols += [snv.AC, snv.AN, snv.AF]
                        csv_file.write('\t'.join([str(col) for col in cols]) + '\n')

if __name__ == '__main__':
    main()