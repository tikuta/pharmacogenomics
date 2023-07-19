#!/usr/bin/env python3
import os
import glob
import json

save_dir = './SNVs-with-generic-numbers'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

CODON_TABLE = {
        'TTT': 'Phe', 'TCT': 'Ser', 'TAT': 'Tyr', 'TGT': 'Cys',
        'TTC': 'Phe', 'TCC': 'Ser', 'TAC': 'Tyr', 'TGC': 'Cys',
        'TTA': 'Leu', 'TCA': 'Ser', 'TAA': '*',   'TGA': '*',
        'TTG': 'Leu', 'TCG': 'Ser', 'TAG': '*',   'TGG': 'Trp',

        'CTT': 'Leu', 'CCT': 'Pro', 'CAT': 'His', 'CGT': 'Arg',
        'CTC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg',
        'CTA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg',
        'CTG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',

        'ATT': 'Ile', 'ACT': 'Thr', 'AAT': 'Asn', 'AGT': 'Ser',
        'ATC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser',
        'ATA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg',
        'ATG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',

        'GTT': 'Val', 'GCT': 'Ala', 'GAT': 'Asp', 'GGT': 'Gly',
        'GTC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly',
        'GTA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly',
        'GTG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'
}

class CDS:
    """
    Note that this class is unconscious of chromosomes  
    """
    def __init__(self, fname) -> None:
        self.genomic_ranges = []
        self.cds_ranges = []
        self.seqs = []
        with open(fname) as f:
            for i, r in enumerate(json.load(f)):
                cols = r['id'].split(':')
                self.chr = cols[2]
                self.genomic_ranges.append([int(cols[3]), int(cols[4])])
                
                seq = r['seq']
                self.seqs.append(r['seq'])

                prev_cds_end = self.cds_ranges[i - 1][1] if i > 0 else 0
                cds_start = prev_cds_end + 1 # 1-origin, inclusive
                cds_end = cds_start + len(seq) - 1 # inclusive
                self.cds_ranges.append([cds_start, cds_end])

    def genomic_to_cds_coord(self, genomic_coord) -> int:
        for genomic_range, cds_range in zip(self.genomic_ranges, self.cds_ranges):
            if genomic_range[0] <= genomic_coord <= genomic_range[1]:
                delta = genomic_coord - genomic_range[0]
                cds_coord = cds_range[0] + delta
                return cds_coord
        raise Exception("Genomic coordinate {} is not included in the CDS ranges ({}).".format(genomic_coord, '/'.join(['-'.join([str(v) for v in c]) for c in self.genomic_ranges])))
    
    def base_at_cds_coord(self, cds_coord) -> str:
        for cds_range, seq in zip(self.cds_ranges, self.seqs):
            if cds_range[0] <= cds_coord <= cds_range[1]:
                idx = cds_coord - cds_range[0]
                base = seq[idx]
                return base
        raise Exception("CDS coordinate {} is not included in the CDS ({}-{}).".format(cds_coord, self.cds_ranges[0][0], self.cds_ranges[-1][1]))
    
    def __str__(self) -> str:
        s = ''
        for genomic_range, cds_range in zip(self.genomic_ranges, self.cds_ranges):
            s += "{} ({})".format(genomic_range[0], cds_range[0])
            s += " - "
            s += "{} ({})".format(genomic_range[1], cds_range[1])
            s += " / "
        return s[:-3]


def main():
    for txt in glob.glob("./SNVs/*.txt"):
        gene = os.path.splitext(os.path.basename(txt))[0]
        w = open(os.path.join(save_dir, gene + ".csv"), 'w')

        cds = CDS("../EnsEMBL/sequence/{}.json".format(gene))
        print(gene, cds)

        generic_nums = {}
        gpcrdb_nums = {}
        with open("../GPCRdb/align/{}.aln".format(gene)) as f:
            for l in f.readlines():
                if l.startswith("#") or len(l) == 0:
                    continue

                _, gpcrdb_num, generic_num, _, ensembl_num = l.strip().split('\t')
                if ensembl_num != '-':
                    generic_nums[int(ensembl_num)] = generic_num
                    gpcrdb_nums[int(ensembl_num)] = gpcrdb_num

        with open(txt) as f:
            # header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "AC", "AN", "AF", "AC_XX", "AN_XX", "AF_XX", "AC_XY", "AN_XY", "AF_XY"]
            for l in f.readlines():
                if len(l) == 0:
                    continue
                if l.startswith("#"):
                    header = l.strip().split('\t')
                    new_header = header[:4] + ["REF_CODON", "REF_AA", "GENERIC_NUM", "ALT_AA", "ALT_CODON"] + header[4:]
                    w.write('\t'.join(new_header) + '\n')
                    continue
                cols = l.strip().split('\t')

                pos = int(cols[1])
                ref = cols[3]
                alt = cols[4]

                cds_coord = cds.genomic_to_cds_coord(pos)
                if ref != cds.base_at_cds_coord(cds_coord):
                    print(pos, cds_coord, cds.base_at_cds_coord(cds_coord))
                    print(ref, alt, cds_coord % 3)
                aa_num, frame = int(cds_coord / 3) + 1, cds_coord % 3
                #generic_num = generic_nums[aa_num]
                
                if frame == 1: # SNV at the first base
                    coords = [cds_coord, cds_coord + 1, cds_coord + 2]
                elif frame == 2: # SNV at the second base
                    coords = [cds_coord - 1, cds_coord, cds_coord + 1]
                else: # SNV at the third base
                    coords = [cds_coord - 2, cds_coord - 1, cds_coord]

                bases_ref = [cds.base_at_cds_coord(coord) for coord in coords]
                codon_ref = ''.join(bases_ref).upper()
                aa_ref = CODON_TABLE[codon_ref]
                
                bases_alt = bases_ref
                if frame == 1:
                    bases_alt[0] = alt
                elif frame == 2:
                    bases_alt[1] = alt
                else:
                    bases_alt[2] = alt
                codon_alt = ''.join(bases_alt).upper()
                aa_alt = CODON_TABLE[codon_alt]

                #new_cols = cols[:4] + [codon_ref, aa_ref, generic_num, aa_alt, codon_alt] + cols[4:]
                new_cols = cols[:4] + [codon_ref, aa_ref, aa_alt, codon_alt, str(aa_num), gpcrdb_nums[aa_num]] + cols[4:]
                w.write('\t'.join(new_cols) + '\n')
        w.close()

if __name__ == '__main__':
    main()