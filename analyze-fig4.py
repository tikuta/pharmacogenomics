#!/usr/bin/env python3
import gpcrdb
import matplotlib
matplotlib.rc('pdf', fonttype=42)
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = "Arial"
from utils import VariationType
import numpy as np
import ensembl
import json
from config import AM_THRESHOLD_BENIGN, AM_THRESHOLD_PATHOGENIC

def analyze_jpn_vs_global(filename_A, filename_B, filename_C, filename_D):
    jpn_vars = {}
    global_vars = {}
    for receptor in gpcrdb.get_filtered_receptor_list():
        with open(receptor.ensembl_path) as f:
            display_name = json.load(f)['display_name']
        with open(receptor.japan_cds_csv_path) as f:
            for l in f.readlines():
                try:
                    anno = ensembl.Annotation.from_csv_line(l)                    
                except ensembl.BlankLineError:
                    continue

                if anno.var_type != VariationType.MISSENSE:
                    continue

                if not anno.pathogenicity:
                    continue

                jpn_vars[(display_name, anno.snv.chromosome, anno.snv.position, anno.snv.ref, anno.snv.alt)] = anno

        with open(receptor.global_cds_csv_path) as f:
            for l in f.readlines():
                try:
                    anno = ensembl.Annotation.from_csv_line(l)                    
                except ensembl.BlankLineError:
                    continue

                if anno.var_type != VariationType.MISSENSE:
                    continue

                if not anno.pathogenicity:
                    continue

                global_vars[(display_name, anno.snv.chromosome, anno.snv.position, anno.snv.ref, anno.snv.alt)] = anno

    all_var_keys = set(jpn_vars.keys()) | set(global_vars.keys())
    pathos = {k: jpn_vars[k].pathogenicity if k in jpn_vars.keys() else global_vars[k].pathogenicity for k in all_var_keys}
    
    likely_pathogenic_keys = tuple(k for k, v in pathos.items() if AM_THRESHOLD_PATHOGENIC < v)
    ambigous_keys = tuple(k for k, v in pathos.items() if AM_THRESHOLD_BENIGN <= v <= AM_THRESHOLD_PATHOGENIC)
    likely_benign_keys = tuple(k for k, v in pathos.items() if v < AM_THRESHOLD_BENIGN)

    likely_pathogenic_freqs = {k: (jpn_vars[k].snv.AF if k in jpn_vars.keys() else 0, global_vars[k].snv.AF if k in global_vars.keys() else 0) for k in likely_pathogenic_keys}
    ambigous_freqs = {k: (jpn_vars[k].snv.AF if k in jpn_vars.keys() else 0, global_vars[k].snv.AF if k in global_vars.keys() else 0) for k in ambigous_keys}
    likely_benign_freqs = {k: (jpn_vars[k].snv.AF if k in jpn_vars.keys() else 0, global_vars[k].snv.AF if k in global_vars.keys() else 0) for k in likely_benign_keys}

    # For Fig. 4a
    fig, ax = plt.subplots(1, 1, figsize=(4, 4), dpi=300)

    ax.fill_between([-1, 2], [-0.75, 2.25], [-1.25, 1.75], color='lightgray', lw=0, alpha=0.8, zorder=-1000)
    ax.fill_between([-1, 2], [-0.5, 2.5], [-1.5, 1.5], color='lightgray', lw=0, alpha=0.6, zorder=-900)
    ax.fill_between([-1, 2], [-0.25, 2.75], [-1.75, 1.25], color='lightgray', lw=0, alpha=0.4, zorder=-800)

    # AF <= 0.25 will be plotted in another fig
    ax.fill_between([0, 0.25], [0, 0], [0.25, 0.25], color='white', lw=0, zorder=-600)
    ax.text(0.125, 0.125, "Fig. S4a", va='center', ha='center', color='lightgray', zorder=-500)

    # Plot only AF > 0.25
    frequent_likely_benign_freqs = [freqs for freqs in likely_benign_freqs.values() if freqs[0] > 0.25 or freqs[1] > 0.25]
    frequent_ambigous_freqs = [freqs for freqs in ambigous_freqs.values() if freqs[0] > 0.25 or freqs[1] > 0.25]
    frequent_likely_pathogenic_freqs = [freqs for freqs in likely_pathogenic_freqs.values() if freqs[0] > 0.25 or freqs[1] > 0.25]
    ax.scatter(*np.array(frequent_likely_benign_freqs).T, marker='.', color='tab:gray', label="likely benign", zorder=1)
    ax.scatter(*np.array(frequent_ambigous_freqs).T, marker='.', color='tab:orange', label="ambigous", zorder=10)
    ax.scatter(*np.array(frequent_likely_pathogenic_freqs).T, marker='D', s=10, color='tab:orange', label="likely pathogenic", zorder=100)

    notables = []
    for k in likely_pathogenic_keys + ambigous_keys:
        freqs = (jpn_vars[k].snv.AF if k in jpn_vars.keys() else 0, global_vars[k].snv.AF if k in global_vars.keys() else 0)
        if freqs[0] > 0.25 or freqs[1] > 0.25:
            anno: ensembl.Annotation = jpn_vars.get(k, global_vars[k])
            display_name = k[0]
            superscript = anno.generic_number if anno.generic_number else anno.segment
            substitution = "{}{}$^{{{}}}${}".format(anno.ref_aa, anno.residue_number, superscript, anno.alt_aa)
            notables.append([display_name, substitution, anno.pathogenicity, *freqs])

    with open(filename_A, 'w') as f:
        f.write('\n'.join([','.join([str(v) for v in l]) for l in notables]))

    ax.legend()
    lims = (-0.02, 1.02)
    ticks = [0, 0.25, 0.5, 0.75, 1]
    ax.set_xlim(*lims)
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticks)
    ax.set_yticks(ticks)
    ax.set_yticklabels(ticks)
    ax.set_ylim(*lims)
    ax.set_aspect('equal')
    ax.set_xlabel("Allele Freq. in 54KJPN")
    ax.set_ylabel("Allele Freq. in 1KGP")
    fig.tight_layout()
    fig.savefig(filename_B)
    plt.close(fig)

    # Fig. S4a
    fig, ax = plt.subplots(1, 1, figsize=(6, 4), dpi=300)
    ax.set_facecolor(color='#d4d4d4')

    # Plot only AF <= 0.25
    infrequent_likely_benign_freqs = [freqs for freqs in likely_benign_freqs.values() if freqs[0] <= 0.25 and freqs[1] <= 0.25]
    infrequent_ambigous_freqs = [freqs for freqs in ambigous_freqs.values() if freqs[0] <= 0.25 and freqs[1] <= 0.25]
    infrequent_likely_pathogenic_freqs = [freqs for freqs in likely_pathogenic_freqs.values() if freqs[0] <= 0.25 and freqs[1] <= 0.25]
    ax.scatter(*np.array(infrequent_likely_benign_freqs).T, marker='.', color='tab:gray', label="likely benign", zorder=1)
    ax.scatter(*np.array(infrequent_ambigous_freqs).T, marker='.', color='tab:orange', label="ambigous", zorder=10)
    ax.scatter(*np.array(infrequent_likely_pathogenic_freqs).T, marker='D', s=10, color='tab:orange', label="likely pathogenic", zorder=100)

    notables = []
    for k in likely_pathogenic_keys + ambigous_keys:
        freqs = (jpn_vars[k].snv.AF if k in jpn_vars.keys() else 0, global_vars[k].snv.AF if k in global_vars.keys() else 0)
        if freqs[0] <= 0.25 and freqs[1] <= 0.25:
            if 0.1 < freqs[0] or 0.1 < freqs[1]:
                anno: ensembl.Annotation = jpn_vars.get(k, global_vars[k])
                display_name = k[0]
                superscript = anno.generic_number if anno.generic_number else anno.segment
                substitution = "{}{}$^{{{}}}${}".format(anno.ref_aa, anno.residue_number, superscript, anno.alt_aa)
                notables.append([display_name, substitution, anno.pathogenicity, *freqs])

    with open(filename_C, 'w') as f:
        f.write('\n'.join([','.join([str(v) for v in l]) for l in notables]))

    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    lims = (-0.01, 0.26)
    ticks = [0, 0.05, 0.1, 0.15, 0.2, 0.25]
    ax.set_xlim(*lims)
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticks)
    ax.set_yticks(ticks)
    ax.set_yticklabels(ticks)
    ax.set_ylim(*lims)
    ax.set_aspect('equal')
    ax.set_xlabel("Allele Freq. in 54KJPN")
    ax.set_ylabel("Allele Freq. in 1KGP")
    fig.tight_layout()
    fig.savefig(filename_D)
    plt.close(fig)

if __name__ == '__main__':
    analyze_jpn_vs_global("./figures/4b_notables.csv", "./figures/4a_jpn_vs_global.pdf", "./figures/S4a_notables.csv", "./figures/S4a_jpn_vs_global.pdf")