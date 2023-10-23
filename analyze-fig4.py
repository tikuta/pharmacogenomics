#!/usr/bin/env python3
from misc import CSV_JPN_CDS_FILENAME, CSV_GLOBAL_CDS_FILENAME
import gpcrdb
import matplotlib
matplotlib.use('Agg')
matplotlib.rc('pdf', fonttype=42)
import matplotlib.pyplot as plt
from matplotlib import patches
plt.rcParams['font.family'] = "Arial"
import os
from utils import VariationType, Segment
import numpy as np

def analyze_jpn_vs_global():
    jpn_vars = {}
    global_vars = {}
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
                    mut_name = cols[header.index('Ref_AA')] + cols[header.index('Res_Num')] + cols[header.index('Alt_AA')]
                    key = entry_name + " " + mut_name
                    freq = float(cols[header.index('AF')])
                    jpn_vars[key] = freq

        with open(os.path.join(dpath, CSV_GLOBAL_CDS_FILENAME)) as f:
            for l in f.readlines():
                cols = l.strip('\n').split('\t')
                if l.startswith('#'):
                    header = cols
                    continue
                var_type = cols[header.index('#Var_Type')]
                if var_type == VariationType.MISSENSE.name:
                    mut_name = cols[header.index('Ref_AA')] + cols[header.index('Res_Num')] + cols[header.index('Alt_AA')]
                    key = entry_name + " " + mut_name
                    freq = float(cols[header.index('AF')])
                    global_vars[key] = freq
    
    jpn_high_freq_vars = {k for k, v in jpn_vars.items() if v > 0.5}
    global_high_freq_vars = {k for k, v in global_vars.items() if v > 0.5}
    jpn_all_vars = set(jpn_vars.keys())
    global_all_vars = set(global_vars.keys())

    print("AF > 0.5 (Fig. 4a)")
    print("JPN", len(jpn_high_freq_vars))
    print("Global", len(global_high_freq_vars))
    common = jpn_high_freq_vars & global_high_freq_vars
    print("JPN ^ Global", len(common))
    print("only in JPN", len(jpn_high_freq_vars - common))
    print("only in Global", len(global_high_freq_vars - common))
    print('\n'.join([v.replace("_human", "").upper() + " ({}, {})".format(jpn_vars.get(v, 0), global_vars.get(v, 0)) for v in sorted(list(global_high_freq_vars - common), key=lambda v: (jpn_vars.get(v, 0), global_vars.get(v, 0)))]))

    # For Fig. 4a
    fig, ax = plt.subplots(1, 1, figsize=(5, 4), dpi=300)

    only_in_jpn = np.array([(jpn_vars.get(k, 0), global_vars.get(k, 0)) for k in jpn_high_freq_vars - global_high_freq_vars])
    only_in_global = np.array([(jpn_vars.get(k, 0), global_vars.get(k, 0)) for k in global_high_freq_vars - jpn_high_freq_vars])
    both_in_jpn_and_global = np.array([(jpn_vars.get(k, 0), global_vars.get(k, 0)) for k in common])

    ax.axvspan(0.5, 1.05, color='lightgray', alpha=0.6, lw=0)
    ax.axhspan(0.5, 1.05, color='lightgray', alpha=0.6, lw=0)

    ax.scatter(both_in_jpn_and_global.T[0], both_in_jpn_and_global.T[1], marker='o', ec='darkgray', fc='white', label="AF > 0.5")
    ax.scatter(only_in_jpn.T[0], only_in_jpn.T[1], marker='v', color='tab:orange', label="AF > 0.5  in 54KJPN\nAF $\leq$ 0.5  in 1KGP")
    ax.scatter(only_in_global.T[0], only_in_global.T[1], marker='^', color='tab:gray', label="AF > 0.5  in 1KGP\nAF $\leq$ 0.5  in 54KJPN")

    ax.legend(bbox_to_anchor=(0.25, 0.25), loc='center', fontsize=9, ncol=1, edgecolor='white', borderpad=0, labelspacing=1)

    ax.set_xlim(-0.05, 1.05)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels([0, 0.25, 0.5, 0.75, 1], size=9)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticklabels([0, 0.25, 0.5, 0.75, 1], size=9)
    ax.set_ylim(-0.05, 1.05)
    ax.set_aspect('equal')
    ax.set_xlabel("Allele Freq. in 54KJPN", size=9)
    ax.set_ylabel("Allele Freq. in 1KGP", size=9)
    fig.tight_layout()
    fig.savefig("./figures/4a_jpn_vs_global.pdf")

    notable_vars = [k for k in jpn_high_freq_vars - global_high_freq_vars if global_vars.get(k, 0) < 0.05]
    notable_vars.sort(key=lambda v: jpn_vars[v], reverse=True)

    for v in [(k, jpn_vars[k], global_vars.get(k, 0)) for k in notable_vars]:
        text = v[0].replace("_human", "").upper()
        print(text, v[1], v[2])

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
                    mut_name = cols[header.index('Ref_AA')] + cols[header.index('Res_Num')] + cols[header.index('Alt_AA')]
                    key = entry_name + " " + mut_name
                    freq = float(cols[header.index('AF')])
                    try:
                        patho = float(cols[header.index('pathogenicity')])
                    except:
                        patho = -1
                    jpn_vars[key] = (freq, patho)

    fig, ax = plt.subplots(1, 1, figsize=(5, 4), dpi=300)
    ax.axhspan(0.34, 0.564, color='gainsboro')
    ax.text(1.01, 0.34 / 2, "likely benign", rotation=90, color="darkgray", va='center', ha='left', size=9)
    ax.text(1.01, (0.564 + 0.34) / 2, "ambiguous", rotation=90, color="darkgray", va='center', ha='left', size=9)
    ax.text(1.01, (0.564 + 1) / 2, "likely pathogenic", rotation=90, color="darkgray", va='center', ha='left', size=9)

    mut_names = list(jpn_vars.keys())
    freqs = [jpn_vars[k][0] for k in mut_names]
    pathos = [jpn_vars[k][1] for k in mut_names]
    ax.scatter(freqs, pathos, marker='.', color='tab:orange')
    for mut_name, freq, patho in zip(mut_names, freqs, pathos):
        if freq > 0.25 and patho > 0.5:
            ax.text(freq - 0.005, patho + 0.005, mut_name.replace(" ", "\n").replace("_human", "").upper(), ha='right', va='bottom', size=7)

    ax.set_xlim(0, 1.05)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels([0, 0.25, 0.5, 0.75, 1], size=9)
    ax.set_yticks([0, 0.34, 0.564, 1])
    ax.set_yticklabels([0, 0.34, 0.564, 1], size=9)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Allele Freq. in 54KJPN", size=9)
    ax.set_ylabel("AlphaMissense pathogenicity", size=9)
    fig.tight_layout()
    fig.savefig("./figures/S4a_pathogenicity.pdf")

if __name__ == '__main__':
    #analyze_jpn_vs_global()
    analyze_alpha_missense_pathogenicity()