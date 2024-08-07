VCF_JPN_GENE_FILENAME = "JPN-GENE.vcf"
VCF_JPN_CDS_FILENAME = "JPN-CDS.vcf"
CSV_JPN_CDS_FILENAME = "JPN-CDS.csv"

VCF_GLOBAL_GENE_FILENAME = "GLOBAL-GENE.vcf"
VCF_GLOBAL_CDS_FILENAME = "GLOBAL-CDS.vcf"
CSV_GLOBAL_CDS_FILENAME = "GLOBAL-CDS.csv"

ALIGNMENT_CANDIDATES_FILENAME = "alignment-candidates.txt"

AM_THRESHOLD_PATHOGENIC = 0.564
AM_THRESHOLD_BENIGN = 0.34
AM_FILENAME = "alpha_missense.tsv"

REFERENCE_CHIMPANZEE_GENOME = "Pan_tro_3.0"

BLOCKLIST = (
    # pseudogenes
    'npy6r_human', # https://www.uniprot.org/uniprotkb/Q99463/entry
    'taar3_human', # https://www.uniprot.org/uniprotkb/Q9P1P4/entry
    'agre4_human', # https://www.uniprot.org/uniprotkb/Q86SQ3/entry
    'agrf2_human', # https://www.uniprot.org/uniprotkb/Q8IZF7/entry
    't2r45_human', # https://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000261936;r=HSCHR12_2_CTG2:327495-328424;t=ENST00000571573
    # inconsistent sequences between GPCRdb and Ensembl
    'gp107_human',
    'g137a_human',
    # olfactory receptors
    'o51i1_human', 
    'o51e1_human',
    'or1g1_human',
    'o51j1_human',
    'o51s1_human',
    'o51a2_human',
    'o2t11_human',
    'o51a4_human',
    'o51b6_human',
    'o51t1_human',
    'o51d1_human',
    'o51e2_human',
    'o51b5_human',
    'o51b2_human',
    'o51l1_human',
    'o51g2_human',
    'o51b4_human',
    'o51a7_human',
    'o51q1_human',
    'o51f2_human',
    'or1a1_human',
    'o51g1_human',
    'o51m1_human',
    'o51i2_human',
    'o51f1_human',
    'o51v1_human'
)

CHIMPANZEE_BLOCK_LIST = (
    # These canonical transcripts did not start with ATG.
    'ENSPTRG00000045080',
    'ENSPTRG00000019404', 
    'ENSPTRG00000019043', 
    'ENSPTRG00000052258', 
    'ENSPTRG00000049073', 
    'ENSPTRG00000009607', 
    'ENSPTRG00000045634', 
    'ENSPTRG00000009279',
    # Broken proteins
    'ENSPTRG00000048811'
)

STANDARD_CODES = {
    'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
    'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
    'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
    'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',

    'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
    'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
    'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
    'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',

    'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
    'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
    'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
    'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',

    'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
    'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
    'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
    'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'
}