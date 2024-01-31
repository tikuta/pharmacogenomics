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

def _is_N_glycosylation_motif(triplet):
    assert(len(triplet) == 3)

    # Check if the variation is related with [N]-X-S/T motif, where X is not Proline.
    if triplet[0] == 'N' and triplet[1] != 'P' and triplet[2] in 'ST':
        return True
    return False

def analyze_terminal_regions():
    n_glyco_gain, n_glyco_loss = [], []
    o_glyco_gain, o_glyco_loss = [], []
    phospho_gain, phospho_loss = [], []
    for receptor in gpcrdb.get_filtered_receptor_list("receptors.json"):
        with open(receptor.ensembl_path) as f:
            display_name = json.load(f)['display_name']
        with open(receptor.alignment_path) as f:
            seq = ''.join([l.split(',')[0][0] for l in f if not l.startswith('#')])
        with open(receptor.japan_cds_csv_path) as f:
            for l in f:
                try:
                    anno = ensembl.Annotation.from_csv_line(l)  
                except ensembl.BlankLineError:
                    continue
                
                if anno.var_type != VariationType.MISSENSE:
                    continue
                if anno.snv.AF <= 0.5:
                    continue
                assert(seq[anno.residue_number - 1] == anno.ref_aa)
                
                # Check glycosylation
                if anno.segment == Segment.Nterm:    
                    # N-glycosylation
                    if 0 <= anno.residue_number - 3 and anno.residue_number <= len(seq):
                        ref_leading_triplet = seq[anno.residue_number - 3:anno.residue_number]
                        alt_leading_triplet = ref_leading_triplet[0] + ref_leading_triplet[1]+ anno.alt_aa

                        if _is_N_glycosylation_motif(ref_leading_triplet) != _is_N_glycosylation_motif(alt_leading_triplet):
                            if _is_N_glycosylation_motif(ref_leading_triplet):
                                # Loss
                                n_glyco_loss.append({"display_name": display_name, "annotation": anno, 
                                                    "ref_motif": ref_leading_triplet, "alt_motif": alt_leading_triplet})
                            else:
                                # Gain
                                n_glyco_gain.append({"display_name": display_name, "annotation": anno,
                                                    "ref_motif": ref_leading_triplet, "alt_motif": alt_leading_triplet})

                    if 0 <= anno.residue_number - 2 and anno.residue_number + 1 <= len(seq):
                        ref_midmost_triplet = seq[anno.residue_number - 2:anno.residue_number + 1]
                        alt_midmost_triplet = ref_midmost_triplet[0] + anno.alt_aa + ref_midmost_triplet[2]

                        if _is_N_glycosylation_motif(ref_midmost_triplet) != _is_N_glycosylation_motif(alt_midmost_triplet):
                            if _is_N_glycosylation_motif(ref_midmost_triplet):
                                # Loss
                                n_glyco_loss.append({"display_name": display_name, "annotation": anno,
                                                    "ref_motif": ref_midmost_triplet, "alt_motif": alt_midmost_triplet})
                            else:
                                # Gain
                                n_glyco_gain.append({"display_name": display_name, "annotation": anno,
                                                    "ref_motif": ref_midmost_triplet, "alt_motif": alt_midmost_triplet})

                    if 0 <= anno.residue_number - 1 and anno.residue_number + 2 <= len(seq):
                        ref_tailing_triplet = seq[anno.residue_number - 1:anno.residue_number + 2]
                        alt_tailing_triplet = anno.alt_aa + ref_tailing_triplet[1] + ref_tailing_triplet[2]

                        if _is_N_glycosylation_motif(ref_tailing_triplet) != _is_N_glycosylation_motif(alt_tailing_triplet):
                            if _is_N_glycosylation_motif(ref_tailing_triplet):
                                # Loss
                                n_glyco_loss.append({"display_name": display_name, "annotation": anno,
                                                    "ref_motif": ref_tailing_triplet, "alt_motif": alt_tailing_triplet})
                            else:
                                # Gain
                                n_glyco_gain.append({"display_name": display_name, "annotation": anno,
                                                    "ref_motif": ref_tailing_triplet, "alt_motif": alt_tailing_triplet})
                    # O-glycosylation
                    if anno.ref_aa in 'ST' and anno.alt_aa not in 'ST':
                        # Loss
                        o_glyco_loss.append({"display_name": display_name, "annotation": anno})
                    elif anno.ref_aa not in 'ST' and anno.alt_aa in 'ST':
                        # Gain
                        o_glyco_gain.append({"display_name": display_name, "annotation": anno})
                # Check phosphorylation
                if anno.segment == Segment.Cterm:
                    if anno.ref_aa in 'ST' and anno.alt_aa not in 'ST':
                        # Loss
                        phospho_loss.append({"display_name": display_name, "annotation": anno})
                    elif anno.ref_aa not in 'ST' and anno.alt_aa in 'ST':
                        # Gain
                        phospho_gain.append({"display_name": display_name, "annotation": anno})
    n_glyco_gain.sort(key=lambda d: d['annotation'].snv.AF)
    n_glyco_loss.sort(key=lambda d: d['annotation'].snv.AF)
    o_glyco_gain.sort(key=lambda d: d['annotation'].snv.AF)
    o_glyco_loss.sort(key=lambda d: d['annotation'].snv.AF)
    phospho_gain.sort(key=lambda d: d['annotation'].snv.AF)
    phospho_loss.sort(key=lambda d: d['annotation'].snv.AF)

    fig, axes = plt.subplots(3, 1, figsize=(4, 6), dpi=300, sharex=True, sharey=True)
    for ax in axes:
        ax.set_facecolor('whitesmoke')
        ax.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1])
        ax.set_yticklabels([0.5, 0.6, 0.7, 0.8, 0.9, 1])
        ax.set_ylabel("Allele Freq.")

    ax = axes[0]
    n_glyco = [[d['annotation'].snv.AF for d in n_glyco_gain], [d['annotation'].snv.AF for d in n_glyco_loss]]
    ax.hist(n_glyco, bins=10, range=(0.5, 1), stacked=True, label=["Gain", "Loss"], orientation='horizontal',
            color=[Segment.Nterm.color, 'white'], edgecolor=Segment.Nterm.color)
    ax.set_title("N-glycosylation (Asn-X-Ser/Thr)") 
    for glyco in n_glyco_gain + n_glyco_loss:
        anno = glyco['annotation']
        label = "{} {}{}{} ({}$\\rightarrow${})".format(glyco['display_name'], anno.ref_aa, anno.residue_number, anno.alt_aa,
                                              glyco['ref_motif'], glyco['alt_motif'])
        ranges = [[0.8, 0.85], [0.85, 0.9], [0.9, 0.95], [0.95, 1]]
        for r in ranges:
            if r[0] < anno.snv.AF <= r[1]:
                ax.text(1, (r[0] + r[1]) / 2, " " + label, ha='left', va='center', size=9)
    ax.legend(loc='lower right', ncol=2)

    ax = axes[1]
    o_glyco = [[d['annotation'].snv.AF for d in o_glyco_gain], [d['annotation'].snv.AF for d in o_glyco_loss]]
    ax.hist(o_glyco, bins=10, range=(0.5, 1), stacked=True, label=["Gain", "Loss"], orientation='horizontal',
            color=[Segment.Nterm.color, 'white'], edgecolor=Segment.Nterm.color)
    ax.set_title("O-glycosylation (Ser/Thr)")
    
    ax = axes[2]
    phospho = [[d['annotation'].snv.AF for d in phospho_gain], [d['annotation'].snv.AF for d in phospho_loss]]
    ax.hist(phospho, bins=10, range=(0.5, 1), stacked=True, label=["Gain", "Loss"], orientation='horizontal',
            color=[Segment.Cterm.color, 'white'], edgecolor=Segment.Cterm.color)
    ax.set_title("Phosphorylation (Ser/Thr)")
    ax.legend(loc='lower right', ncol=2)
    ax.set_xticks(range(5))
    ax.set_xticklabels(range(5))
    ax.set_xlabel("Number of Variants")

    fig.tight_layout()    
    fig.savefig("./figures/2c_term.pdf")

if __name__ == '__main__':
    analyze_high_allele_freq_vars()
    analyze_terminal_regions()