#!/usr/bin/env python3
from misc import CSV_JPN_CDS_FILENAME
import gpcrdb
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = "Arial"
import os
from utils import VariationType, Segment
import json

def analyze_high_allele_freq_vars():
    vars = []
    for receptor in gpcrdb.get_filtered_receptor_list("receptors.json"):
        entry_name = receptor['entry_name']
        receptor_class = receptor['receptor_class']
        dpath = os.path.join(receptor_class, entry_name)
        with open(os.path.join(dpath, 'ensembl.json')) as f:
            display_name = json.load(f)['display_name']
        with open(os.path.join(dpath, CSV_JPN_CDS_FILENAME)) as f:
            for l in f.readlines():
                cols = l.strip('\n').split('\t')
                if l.startswith('#'):
                    header = cols
                    continue
                var_type = cols[header.index('#Var_Type')]
                seg = Segment.value_of(cols[header.index('Seg')])
                freq = float(cols[header.index('AF')])
                aa_ref = cols[header.index('Ref_AA')]
                res_num = cols[header.index('Res_Num')]
                aa_alt = cols[header.index('Alt_AA')]
                gen_num = cols[header.index('Generic_Num')]
                if seg == 'None':
                    print(dpath, l)
                if var_type == VariationType.MISSENSE.name and freq > 0.5:
                    label = "{} {}{}{} ({})".format(display_name, aa_ref, res_num, aa_alt, gen_num if gen_num != 'None' else seg)
                    classA = True if receptor_class == 'Class A (Rhodopsin)' else False
                    var = [freq, seg, label, classA]
                    vars.append(var)
    vars.sort(key=lambda v: (v[0], -1 * v[1].index))

    # Fig. 2a
    fig, ax = plt.subplots(1, 1, figsize=(5, 1.5), dpi=300)
    ax.set_facecolor('whitesmoke')

    left = 0
    for seg in Segment:
        width = sum([True if var[1] == seg else False for var in vars])
        ax.barh(0, width, height=0.7, left=left, color=seg.color, edgecolor='k', lw=0.2)
        
        label = seg.value
        if label.startswith('TM') or label.startswith('H8') or label.endswith('-term'):
            if label in ('TM1', 'N-term', 'H8', 'C-term'):
                ax.text(left + width / 2, 0.4, label, ha='center', va='bottom', size=6)
            else:
                ax.text(left + width / 2, 0.4, label[-1], ha='center', va='bottom', size=6)
        left += width

    ax.set_ylim(-0.5, 0.7)
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.set_xlabel("Number of Missense Variants (AF > 0.5)")
    fig.tight_layout()
    fig.savefig("./figures/2a_high_allele_freq_vars.pdf")
    plt.close(fig)

    # Fig. 2b
    fig, ax = plt.subplots(1, 1, figsize=(4, 12), dpi=300)
    ax.set_facecolor('whitesmoke')

    for i, var in enumerate(vars):
        ax.barh(i, var[0], color=var[1].color)
        freq = var[0]
        if freq == 1:
            text = "1"
        else:
            text = "> {:.2f}".format(int(freq * 100) / 100)
        ax.text(var[0] + 0.01, i, text, size=5, ha='left', va='center')
    
    A_ticks, A_ticklabels = [], []
    non_A_ticks, non_A_ticklabels = [], []
    for i, v in enumerate(vars):
        if v[3]:
            A_ticks.append(i)
            A_ticklabels.append(v[2])
        else:
            non_A_ticks.append(i)
            non_A_ticklabels.append(v[2])
    ax.set_yticks(A_ticks)
    ax.set_yticklabels(A_ticklabels, size=7)
    ax.set_yticks(non_A_ticks, minor=True)
    ax.set_yticklabels(non_A_ticklabels, size=7, minor=True, color='tab:gray')
    ax.set_xlabel("Allele Freq.")
    ax.set_ylim(-0.8, len(vars) - 0.3)
    ax.set_xticks([0.5, 0.6, 0.8, 1.0])
    ax.set_xticklabels([0.5, 0.6, 0.8, 1.0])
    ax.set_xlim(left=0.5, right=1.1)
    fig.tight_layout()
    fig.savefig("./figures/2b_high_allele_freq_vars.pdf")

if __name__ == '__main__':
    analyze_high_allele_freq_vars()