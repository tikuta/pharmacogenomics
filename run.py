#!/usr/bin/env python3

import os
import gpcrdb
import uniprot
import ensembl
import vcf
from utils import Region, Gene, Variation

def main():
    vcfs = [
        "./data/tommo-38kjpn-20220630-GRCh38-af-autosome.vcf.gz", 
        "./data/tommo-38kjpn-20220929-GRCh38-af-chrX_PAR2.vcf.gz",
        "./data/tommo-38kjpn-20220929-GRCh38-af-chrX_PAR3.vcf.gz"
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
        vcf.filter_vcf(vcfs, [gene_region], os.path.join(dpath, "38KJPN-gene.vcf"))

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

            vcf.filter_vcf(vcfs, coding_regions, os.path.join(dpath, "38KJPN-CDS.vcf"))

            ref_seq = gene.translate()
            assert(ref_seq[0] == 'M')
            matched = gpcrdb.match(ref_seq, generic_numbers, receptor_class, os.path.join(dpath, "match.json"), force=True, alignment_for_human=os.path.join(dpath, "match.txt"))

            w = open(os.path.join(dpath, "38KJPN-CDS.csv"), 'w')
            header = ["Var_Type", "Chr", "Pos", "Res_Num"]
            header += ["Ref_Base", "Ref_Codon", "Ref_AA"]
            header += ["Alt_Base", "Alt_Codon", "Alt_AA"]
            header += ["Seg", "Generic_Num"]
            header += ["AC", "AN", "AF"]
            header += ["AC_XX", "AN_XX", "AF_XX"]
            header += ["AC_XY", "AN_XY", "AF_XY"]
            w.write('#' + '\t'.join(header) + '\n')

            with open(os.path.join(dpath, "38KJPN-CDS.vcf")) as f:
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
                            print(anno.res_num, snv.position, snv.ref, snv.alt)
                        assert(res['ensembl_sequence_number'] == anno.res_num)
                        seg = res['segment']
                        generic_num = res['generic_number']
                        cols = [anno.var_type, snv.chromosome, snv.position, anno.res_num]
                        cols += [snv.ref, anno.ref_codon, anno.ref_aa]
                        cols += [snv.alt, anno.alt_codon, anno.alt_aa]
                        cols += [seg, generic_num]
                        cols += [snv.AC, snv.AN, snv.AF]
                        cols += [snv.AC_XX, snv.AN_XX, snv.AF_XX]
                        cols += [snv.AC_XY, snv.AN_XY, snv.AF_XY]
                        w.write('\t'.join([str(col) for col in cols]) + '\n')            

if __name__ == '__main__':
    main()