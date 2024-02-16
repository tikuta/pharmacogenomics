#!/usr/bin/env python3
import gpcrdb
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = "Arial"
import os
from utils import VariationType
import numpy as np
import pandas as pd
import ensembl
import json
from config import AM_THRESHOLD_BENIGN, AM_THRESHOLD_PATHOGENIC

def analyze_jpn_vs_global():
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
    
    likely_pathogenic_keys = (k for k, v in pathos.items() if AM_THRESHOLD_PATHOGENIC < v)
    ambigous_keys = (k for k, v in pathos.items() if AM_THRESHOLD_BENIGN <= v <= AM_THRESHOLD_PATHOGENIC)
    likely_benign_keys = (k for k, v in pathos.items() if v < AM_THRESHOLD_BENIGN)

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
    fig.savefig("./figures/4a_jpn_vs_global.pdf")
    plt.close(fig)

    # Fig. S4a
    fig, ax = plt.subplots(1, 1, figsize=(4, 4), dpi=300)
    ax.set_facecolor(color='#d4d4d4')

    # Plot only AF <= 0.25
    infrequent_likely_benign_freqs = [freqs for freqs in likely_benign_freqs.values() if freqs[0] <= 0.25 and freqs[1] <= 0.25]
    infrequent_ambigous_freqs = [freqs for freqs in ambigous_freqs.values() if freqs[0] <= 0.25 and freqs[1] <= 0.25]
    infrequent_likely_pathogenic_freqs = [freqs for freqs in likely_pathogenic_freqs.values() if freqs[0] <= 0.25 and freqs[1] <= 0.25]
    ax.scatter(*np.array(infrequent_likely_benign_freqs).T, marker='.', color='tab:gray', label="likely benign", zorder=1)
    ax.scatter(*np.array(infrequent_ambigous_freqs).T, marker='.', color='tab:orange', label="ambigous", zorder=10)
    ax.scatter(*np.array(infrequent_likely_pathogenic_freqs).T, marker='D', s=10, color='tab:orange', label="likely pathogenic", zorder=100)

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
    fig.savefig("./figures/S4a_jpn_vs_global.pdf")
    plt.close(fig)
"""
    # For Fig. 4b
    pathogenic_vars = {
        "GPCR": [], 
        "Family": [],
        "Substitution": [],
        "rsID": [],
        "AF in 54KJPN": [],
        "AF in 1KGP": [],
        "AM pathogenicity": []
    }
    for var in jpn_vars.keys():
        freq_jpn = jpn_vars[var][0]
        freq_global = global_vars.get(var, [0])[0]
        patho = jpn_vars[var][1]
        if patho > AM_THRESHOLD_BENIGN and (freq_jpn > 0.25 or freq_global > 0.25):
            pathogenic_vars['GPCR'].append(var.split(' ')[0].replace("_human", "").upper())
            pathogenic_vars['Family'].append(jpn_vars[var][2])
            pathogenic_vars['Substitution'].append(var.split(' ')[1])
            pathogenic_vars['rsID'].append(jpn_vars[var][3])
            pathogenic_vars['AF in 54KJPN'].append(freq_jpn)
            pathogenic_vars['AF in 1KGP'].append(freq_global)
            pathogenic_vars['AM pathogenicity'].append(patho)
    
    df = pd.DataFrame(pathogenic_vars).sort_values('AF in 54KJPN', ascending=False)
    df.to_csv("./tables/4b_jpn_pathogenic_vars.csv")

    # For Fig. S4a
    lims = (-0.01, 0.26)

    fig, ax = plt.subplots(1, 1, figsize=(5, 5), dpi=300)
    ax.set_facecolor(color='#e5e5e5')

    ax.scatter(likely_benign.T[0], likely_benign.T[1], marker='.', color='tab:gray', label="likely benign")
    ax.scatter(ambigous.T[0], ambigous.T[1], marker='.', color='tab:orange', label="ambigous")
    ax.scatter(likely_pathogenic.T[0], likely_pathogenic.T[1], marker='D', s=10, color='tab:orange', label="likely pathogenic")

    for k, v in jpn_vars.items():
        if v[1] > AM_THRESHOLD_BENIGN:
            x, y = jpn_vars.get(k, [0])[0], global_vars.get(k, [0])[0]
            if not (0 <= x < 0.1 and 0 <= y < 0.1) and not (0.25 < x or 0.25 < y):
                ax.text(x, y - 0.005, k.replace("_human ", "\n").upper(), color='tab:orange', ha='center', va='top', size=9)
                print(k, x, y, v[0])

    ax.legend(bbox_to_anchor=(0.5, 1.02), loc='lower center', fontsize=9, ncol=3)

    ax.set_xlim(*lims)
    ax.set_xticks([0, 0.05, 0.1, 0.15, 0.2, 0.25])
    ax.set_xticklabels([0, 0.05, 0.1, 0.15, 0.2, 0.25], size=9)
    ax.set_yticks([0, 0.05, 0.1, 0.15, 0.2, 0.25])
    ax.set_yticklabels([0, 0.05, 0.1, 0.15, 0.2, 0.25], size=9)
    ax.set_ylim(*lims)
    ax.set_aspect('equal')
    ax.set_xlabel("Allele Freq. in 54KJPN", size=9)
    ax.set_ylabel("Allele Freq. in 1KGP", size=9)
    fig.tight_layout()
    fig.savefig("./figures/S4a_jpn_vs_global.pdf")

    # For Fig. S4b
    pathogenic_vars = {
        "GPCR": [], 
        "Family": [],
        "Substitution": [],
        "rsID": [],
        "AF in 54KJPN": [],
        "AF in 1KGP": [],
        "AM pathogenicity": []
    }
    for var in jpn_vars.keys():
        freq_jpn = jpn_vars[var][0]
        freq_global = global_vars.get(var, [0])[0]
        patho = jpn_vars[var][1]
        if patho > AM_THRESHOLD_BENIGN and (freq_jpn <= 0.25 and freq_global <= 0.25):
            pathogenic_vars['GPCR'].append(var.split(' ')[0].replace("_human", "").upper())
            pathogenic_vars['Family'].append(jpn_vars[var][2])
            pathogenic_vars['Substitution'].append(var.split(' ')[1])
            pathogenic_vars['rsID'].append(jpn_vars[var][3])
            pathogenic_vars['AF in 54KJPN'].append(freq_jpn)
            pathogenic_vars['AF in 1KGP'].append(freq_global)
            pathogenic_vars['AM pathogenicity'].append(patho)
    
    df = pd.DataFrame(pathogenic_vars).sort_values('AF in 54KJPN', ascending=False)
    df.to_csv("./figures/S4b_jpn_pathogenic_vars.csv")

    # For Table S1
    pathogenic_vars = {
        "GPCR": [], 
        "Family": [],
        "Substitution": [],
        "rsID": [],
        "AF in 54KJPN": [],
        "AF in 1KGP": [],
        "AM pathogenicity": []
    }
    for var in jpn_vars.keys():
        freq_jpn = jpn_vars[var][0]
        freq_global = global_vars.get(var, [0])[0]
        patho = jpn_vars[var][1]
        if patho > AM_THRESHOLD_PATHOGENIC:
            pathogenic_vars['GPCR'].append(var.split(' ')[0].replace("_human", "").upper())
            pathogenic_vars['Family'].append(jpn_vars[var][2])
            pathogenic_vars['Substitution'].append(var.split(' ')[1])
            pathogenic_vars['rsID'].append(jpn_vars[var][3])
            pathogenic_vars['AF in 54KJPN'].append(freq_jpn)
            pathogenic_vars['AF in 1KGP'].append(freq_global)
            pathogenic_vars['AM pathogenicity'].append(patho)
    
    df = pd.DataFrame(pathogenic_vars).sort_values('AM pathogenicity', ascending=False)
    df.to_csv("./tables/S1_jpn_pathogenic_vars.csv")
"""
def analyze_alpha_missense_pathogenicity():
    jpn_vars = {}
    for receptor in gpcrdb.get_filtered_receptor_list("receptors.json"):
        entry_name = receptor['entry_name']
        receptor_class = receptor['receptor_class']
        dpath = os.path.join(receptor_class, entry_name)

        with open(os.path.join(dpath, CSV_JPN_CDS_FILENAME)) as f:
            for l in f.readlines():
                cols = l.strip('\n').split('\t')
                if l.startswith('#'):
                    header = cols
                    continue
                var_type = cols[header.index('#Var_Type')]
                if var_type == VariationType.MISSENSE.name:
                    gen_num = cols[header.index('Generic_Num')]
                    if gen_num == 'None':
                        gen_num = cols[header.index('Seg')]
                    mut_name = "{}{}$^{{{}}}${}".format(cols[header.index('Ref_AA')], cols[header.index('Res_Num')], gen_num, cols[header.index('Alt_AA')])
                    key = entry_name + " " + mut_name
                    freq = float(cols[header.index('AF')])
                    try:
                        patho = float(cols[header.index('pathogenicity')])
                    except:
                        patho = -1
                    jpn_vars[key] = (freq, patho)

    fig, ax = plt.subplots(1, 1, figsize=(5, 4), dpi=300)
    ax.axhspan(AM_THRESHOLD_BENIGN, AM_THRESHOLD_PATHOGENIC, color='gainsboro')
    ax.text(1.01, AM_THRESHOLD_BENIGN / 2, "likely benign", rotation=90, color="darkgray", va='center', ha='left', size=9)
    ax.text(1.01, (AM_THRESHOLD_PATHOGENIC + AM_THRESHOLD_BENIGN) / 2, "ambiguous", rotation=90, color="darkgray", va='center', ha='left', size=9)
    ax.text(1.01, (AM_THRESHOLD_PATHOGENIC + 1) / 2, "likely pathogenic", rotation=90, color="darkgray", va='center', ha='left', size=9)

    mut_names = list(jpn_vars.keys())
    freqs = [jpn_vars[k][0] for k in mut_names]
    pathos = [jpn_vars[k][1] for k in mut_names]
    ax.scatter(freqs, pathos, marker='.', color='tab:gray')
    for mut_name, freq, patho in zip(mut_names, freqs, pathos):
        if freq > 0.25 and patho > 0.5:
            ax.text(freq - 0.005, patho + 0.005, mut_name.replace(" ", "\n").replace("_human", "").upper(), ha='right', va='bottom', size=7)

    ax.set_xlim(0, 1.05)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels([0, 0.25, 0.5, 0.75, 1], size=9)
    ax.set_yticks([0, AM_THRESHOLD_BENIGN, AM_THRESHOLD_PATHOGENIC, 1])
    ax.set_yticklabels([0, AM_THRESHOLD_BENIGN, AM_THRESHOLD_PATHOGENIC, 1], size=9)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Allele Freq. in 54KJPN", size=9)
    ax.set_ylabel("AlphaMissense pathogenicity", size=9)
    fig.tight_layout()
    fig.savefig("./figures/S4b_pathogenicity.pdf")

if __name__ == '__main__':
    analyze_jpn_vs_global()
    # analyze_alpha_missense_pathogenicity()