#!/usr/bin/env python3
import gpcrdb
import matplotlib
matplotlib.rc('pdf', fonttype=42)
matplotlib.use('Agg')
from matplotlib import patches
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = "Arial"
from utils import VariationType
import numpy as np
import ensembl
import json
from config import AM_THRESHOLD_BENIGN, AM_THRESHOLD_PATHOGENIC

def analyze_jpn_vs_global(filename_A, filename_B, filename_C, filename_D, filename_E, filename_F):
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
    fig, ax = plt.subplots(1, 1, figsize=(5, 5), dpi=300)

    ax.fill_between([-1, 2], [-0.75, 2.25], [-1.25, 1.75], color='lightgray', lw=0, alpha=0.8, zorder=-1000)
    ax.fill_between([-1, 2], [-0.5, 2.5], [-1.5, 1.5], color='lightgray', lw=0, alpha=0.6, zorder=-900)
    ax.fill_between([-1, 2], [-0.25, 2.75], [-1.75, 1.25], color='lightgray', lw=0, alpha=0.4, zorder=-800)

    # Plot only AF > 0.25
    frequent_likely_benign_freqs = [freqs for freqs in likely_benign_freqs.values() if freqs[0] > 0.25 or freqs[1] > 0.25]
    frequent_ambigous_freqs = [freqs for freqs in ambigous_freqs.values() if freqs[0] > 0.25 or freqs[1] > 0.25]
    frequent_likely_pathogenic_freqs = [freqs for freqs in likely_pathogenic_freqs.values() if freqs[0] > 0.25 or freqs[1] > 0.25]
    ax.scatter(*np.array(frequent_likely_benign_freqs).T, marker='.', color='tab:gray', label="likely benign", zorder=1)
    ax.scatter(*np.array(frequent_ambigous_freqs).T, marker='.', color='tab:orange', label="ambigous", zorder=10)
    ax.scatter(*np.array(frequent_likely_pathogenic_freqs).T, marker='D', s=10, color='tab:orange', label="likely pathogenic", zorder=100)

    # AF <= 0.25 will be plotted in another fig
    num_infrequent = len(all_var_keys) - len(frequent_ambigous_freqs) - len(frequent_likely_benign_freqs) - len(frequent_likely_pathogenic_freqs)
    ax.fill_between([0, 0.25], [0, 0], [0.25, 0.25], color='white', lw=0, zorder=-600)
    ax.text(0.125, 0.125, "Fig. S4a\n{:,} vars.".format(num_infrequent), va='center', ha='center', color='tab:gray', zorder=-500)
    ax.text(0, 1, "Total {:,} vars. (incl. AF$\leq$0.25)".format(len(all_var_keys)), va='top', ha='left', color='tab:gray', zorder=-500)

    # Japanese-specific variants
    jpn_freq_threshold, global_freq_threshold = 0.25, 0.015
    x, y = jpn_freq_threshold, -global_freq_threshold
    w, h = 1 - jpn_freq_threshold + 0.015, global_freq_threshold * 2
    r = patches.Rectangle(xy=(x, y), width=w, height=h, ec='tab:gray', fill=False, zorder=-1, lw=0.6)
    ax.add_patch(r)
    ax.text(1, global_freq_threshold, "Table S1", color='tab:gray', ha='right', va='bottom')

    lines = []
    for k in all_var_keys:
        jpn_var, global_var = jpn_vars.get(k, None), global_vars.get(k, None)
        jpn_freq, global_freq = jpn_var.snv.AF if jpn_var else 0, global_var.snv.AF if global_var else 0
        if jpn_freq > jpn_freq_threshold and global_freq < global_freq_threshold:
            display_name = k[0]
            superscript = jpn_var.generic_number if jpn_var.generic_number else jpn_var.segment
            substitution = "{}{}^{{{}}}{}".format(jpn_var.ref_aa, jpn_var.residue_number, superscript, jpn_var.alt_aa)
            cols = [display_name, substitution, jpn_var.snv.rsid, jpn_freq, global_freq]
            lines.append(cols)
    lines.sort(key=lambda l: float(l[3]), reverse=True)
    with open(filename_E, 'w') as f:
        f.write(",".join(["#GPCR", "Substitution", "rsID", "AF_jpn", "AF_global"]) + "\n")
        f.write('\n'.join([",".join([str(v) for v in l]) for l in lines]))

    notables = []
    for k in likely_pathogenic_keys + ambigous_keys:
        freqs = (jpn_vars[k].snv.AF if k in jpn_vars.keys() else 0, global_vars[k].snv.AF if k in global_vars.keys() else 0)
        if freqs[0] > 0.25 or freqs[1] > 0.25:
            anno: ensembl.Annotation = jpn_vars.get(k, global_vars[k])
            display_name = k[0]
            superscript = anno.generic_number if anno.generic_number else anno.segment
            substitution = "{}{}$^{{{}}}${}".format(anno.ref_aa, anno.residue_number, superscript, anno.alt_aa)
            notables.append([display_name, substitution, anno.snv.rsid, anno.pathogenicity, *freqs])
    notables.sort(key=lambda n: n[5], reverse=True)
    with open(filename_A, 'w') as f:
        f.write('\n'.join([','.join([str(v) for v in l]) for l in notables]))

    ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1), ncol=3)
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
    fig, ax = plt.subplots(1, 1, figsize=(5, 5), dpi=300)
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
                notables.append([display_name, substitution, anno.snv.rsid, anno.pathogenicity, *freqs])
    notables.sort(key=lambda n: n[5], reverse=True)
    with open(filename_C, 'w') as f:
        f.write('\n'.join([','.join([str(v) for v in l]) for l in notables]))

    ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1), ncol=3)
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

    # Fig. S4c
    fig, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=300)

    ax.axhspan(AM_THRESHOLD_BENIGN, AM_THRESHOLD_PATHOGENIC, color='lightgray')

    xys = [(anno.snv.AF, anno.pathogenicity) for anno in jpn_vars.values()]
    ax.scatter(*np.array(xys).T, marker='.', color='tab:gray')

    benign_text = "likely benign\n({:,})".format(len(likely_benign_keys))
    ax.text(1.07, AM_THRESHOLD_BENIGN / 2, benign_text, va='center', ha='center', rotation=90, color='tab:gray', size=8)
    ambigous_text = "ambigous\n({:,})".format(len(ambigous_keys))
    ax.text(1.07, (AM_THRESHOLD_BENIGN + AM_THRESHOLD_PATHOGENIC) / 2, ambigous_text, va='center', ha='center', rotation=90, color='tab:gray', size=8)
    pathogenic_text = "likely pathogenic\n({:,} vars.)".format(len(likely_pathogenic_keys))
    ax.text(1.07, (AM_THRESHOLD_PATHOGENIC + 1) / 2, pathogenic_text, va='center', ha='center', rotation=90, color='tab:gray', size=8)

    ax.set_xlabel("Allel Freq. in 54KJPN")
    ax.set_xlim(-0.02, 1.12)
    ax.set_ylabel("Pathogenicity")
    ax.set_ylim(-0.02, 1.02)
    ax.set_yticks([0, AM_THRESHOLD_BENIGN, AM_THRESHOLD_PATHOGENIC, 1])
    ax.set_yticklabels([0, AM_THRESHOLD_BENIGN, AM_THRESHOLD_PATHOGENIC, 1])

    fig.tight_layout()
    fig.savefig(filename_F)

if __name__ == '__main__':
    analyze_jpn_vs_global("./figures/4a_notables.csv", "./figures/4a_jpn_vs_global.pdf", 
                          "./figures/S4a_notables.csv", "./figures/S4a_jpn_vs_global.png",
                          "./figures/TableS1_japanese_specifics.csv", "./figures/S4c_freq_vs_pathogenicity.png")