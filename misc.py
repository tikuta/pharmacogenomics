import enum
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

class Segment(enum.Enum):
    Nterm = "N-term"
    TM1 = "TM1"
    ICL1 = "ICL1"
    TM2 = "TM2"
    ECL1 = "ECL1"
    TM3 = "TM3"
    ICL2 = "ICL2"
    TM4 = "TM4"
    ECL2 = "ECL2"
    TM5 = "TM5"
    ICL3 = "ICL3"
    TM6 = "TM6"
    ECL3 = "ECL3"
    TM7 = "TM7"
    ICL4 = "ICL4"
    H8 = "H8"
    Cterm = "C-term"
    NONE = "None"

    def __str__(self) -> str:
        return self.value
    
    @property
    def color(self):
        cmap = plt.get_cmap('rainbow', len(Segment) - 1) 
        if self == Segment.NONE:
            return 'tab:gray'
        for i, seg in enumerate(Segment):
            if seg == self:
                return cmap(i)
    
    @classmethod
    def value_of(cls, target):
        for e in cls:
            if e.value == target:
                return e
        raise ValueError

    @classmethod
    def generic_number_of(cls, g_num):
        seg = int(g_num.split('.')[0])
        if 1 <= seg <= 7:
            return cls.value_of("TM" + str(seg))
        elif seg == 8:
            return cls.H8
        elif seg == 12:
            return cls.ICL1
        elif seg == 23:
            return cls.ECL1
        elif seg == 34:
            return cls.ICL2
        elif seg == 45:
            return cls.ECL2
        elif seg == 56:
            return cls.ICL3
        elif seg == 67:
            return cls.ECL3
        raise ValueError

AA2COLOR = {
    "A": "#ffff00",
    "M": "#ffff00",
    "F": "#07b050",
    "Y": "#07b050",
    "W": "#07b050",

    "C": "#bf8f00",
    
    "D": "#ff0000",
    "E": "#ff0000",

    "G": "#ff02ff",

    "H": "#8282D2",

    "I": "#ffff00",
    "L": "#ffff00",
    "V": "#ffff00",

    "K": "#0070c0",
    "R": "#0070c0",

    "N": "#7030a0",
    "Q": "#7030a0",

    "S": "#7030a0",
    "T": "#7030a0",

    "P": "#d603ff",
}


BLOCKLIST = (
    # pseudogenes
    'npy6r_human',
    'taar3_human',
    'agre4_human',
    # ncRNA
    'agrf2_human', # https://www.ncbi.nlm.nih.gov/gene/222611
    # CHR_HSCHR12_2_CTG2
    't2r45_human',
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