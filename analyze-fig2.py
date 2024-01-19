#!/usr/bin/env python3
import gpcrdb
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = "Arial"
import ensembl
from utils import VariationType, Segment
import json

def analyze_high_allele_freq_vars():
    high_frequent_vars = []
    for receptor in gpcrdb.get_filtered_receptor_list("receptors.json"):
        with open(receptor.ensembl_path) as f:
            display_name = json.load(f)['display_name']
        with open(receptor.japan_cds_csv_path) as f:
            for l in f.readlines():
                try:
                    anno = ensembl.Annotation.from_csv_line(l)
                    if anno.segment == Segment.FailedToGuess:
                        if anno.snv.rsid in ('rs376593544', 'rs1249120368'):
                            # Alignment-related issues.
                            # These variations are located between H8 and C-term in CCR2.
                            # (See `Class A (Rhodopsin)/ccr2_human/alignment.csv`)
                            # We reagard these variations as C-term for the following analysis.
                            anno.segment = Segment.Cterm
                        else:
                            raise NotImplementedError
                    elif anno.segment == Segment.NONE:
                        raise NotImplementedError

                    if anno.var_type == VariationType.MISSENSE and anno.snv.AF > 0.5:
                        d = {"display_name": display_name, "annotation": anno, "class A": receptor.receptor_class == 'Class A (Rhodopsin)'}
                        high_frequent_vars.append(d)
                except ensembl.BlankLineError:
                    continue
    high_frequent_vars.sort(key=lambda d: (d['annotation'].snv.AF, -d['annotation'].segment.index, d['class A'], d['annotation'].residue_number))

    # Fig. 2a
    fig, ax = plt.subplots(1, 1, figsize=(4, 1.5), dpi=300)
    ax.set_facecolor('whitesmoke')

    left = 0
    for seg in Segment:
        width = sum([1 if var['annotation'].segment == seg else 0 for var in high_frequent_vars])
        if width == 0:
            continue

        ax.barh(0, width, height=0.7, left=left, color=seg.color, edgecolor='k', lw=0.2)
        
        label = seg.value
        x = left + width / 2
        if label.startswith('TM') or label.startswith('H8') or label.endswith('-term'):
            if label in ('TM1', 'N-term', 'H8', 'C-term'):
                ax.text(x, 0.375, label, ha='center', va='baseline', size=8)
            else:
                ax.text(x, 0.375, label[-1], ha='center', va='baseline', size=8)
        text = str(width) if seg != Segment.Nterm else "{} variants".format(width)
        ax.text(x, -0.375, text, color='black', ha='center', va='top', size=8)
        left += width

    ax.set_xlim(0, left)
    ax.set_ylim(-0.5, 0.5)
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.set_xlabel("Number of Missense Variants (AF > 0.5)")
    fig.tight_layout()
    fig.savefig("./figures/2a_high_allele_freq_vars.pdf")
    plt.close(fig)

    # Fig. 2b
    fig, ax = plt.subplots(1, 1, figsize=(4, 12), dpi=300)
    ax.set_facecolor('whitesmoke')

    for i, var in enumerate(high_frequent_vars):
        anno = var['annotation']
        ax.barh(i, anno.snv.AF, color=anno.segment.color)
        if anno.snv.AF == 1:
            text = "1"
        else:
            text = "> {:.2f}".format(int(anno.snv.AF * 100) / 100)
        ax.text(anno.snv.AF + 0.01, i, text, size=5, ha='left', va='center')

    A_ticks, A_ticklabels = [], []
    non_A_ticks, non_A_ticklabels = [], []
    for i, var in enumerate(high_frequent_vars):
        anno = var['annotation']
        label = "{} {}{}{} ({})".format(var['display_name'], anno.ref_aa, anno.residue_number, anno.alt_aa, 
                                        anno.generic_number if anno.generic_number else anno.segment.value)
        if var['class A']:
            A_ticks.append(i)
            A_ticklabels.append(label)
        else:
            non_A_ticks.append(i)
            non_A_ticklabels.append(label)
    ax.set_yticks(A_ticks)
    ax.set_yticklabels(A_ticklabels, size=7)
    ax.set_yticks(non_A_ticks, minor=True)
    ax.set_yticklabels(non_A_ticklabels, size=7, minor=True, color='tab:gray')
    ax.set_xlabel("Allele Freq.")
    ax.set_ylim(-0.8, len(high_frequent_vars) - 0.3)
    ax.set_xticks([0.5, 0.6, 0.8, 1.0])
    ax.set_xticklabels([0.5, 0.6, 0.8, 1.0])
    ax.set_xlim(left=0.5, right=1.1)
    fig.tight_layout()
    fig.savefig("./figures/2b_high_allele_freq_vars.pdf")

if __name__ == '__main__':
    analyze_high_allele_freq_vars()