#!/usr/bin/env python3
from misc import CSV_JPN_CDS_FILENAME, CSV_GLOBAL_CDS_FILENAME
import gpcrdb
import matplotlib
matplotlib.use('Agg')
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
    
    # For Fig. 4a & S4a
    jpn_high_freq_vars = {k for k, v in jpn_vars.items() if v > 0.5}
    global_high_freq_vars = {k for k, v in global_vars.items() if v > 0.5}
    jpn_all_vars = set(jpn_vars.keys())
    global_all_vars = set(global_vars.keys())

    print("---------------------")
    print("AF > 0.5 (Fig. 4a)")
    print("JPN", len(jpn_high_freq_vars))
    print("Global", len(global_high_freq_vars))
    print("JPN ^ Global", len(jpn_high_freq_vars & global_high_freq_vars))
    print("only in JPN", len(jpn_high_freq_vars - global_high_freq_vars))
    print("only in Global", len(global_high_freq_vars - jpn_high_freq_vars))

    print("---------------------")

    print("All (Fig. S4a)")
    print("JPN", len(jpn_all_vars))
    print("Global", len(global_all_vars))
    print("JPN ^ Global", len(jpn_all_vars & global_all_vars))
    print("only in JPN", len(jpn_all_vars - global_all_vars))
    print("only in Global", len(global_all_vars - jpn_all_vars))
    print("---------------------")

    # For Fig. S4b
    fig, ax = plt.subplots(1, 1, figsize=(6.5, 6), dpi=300)

    only_in_jpn = np.array([(jpn_vars.get(k, 0), global_vars.get(k, 0)) for k in jpn_high_freq_vars - global_high_freq_vars])
    only_in_global = np.array([(jpn_vars.get(k, 0), global_vars.get(k, 0)) for k in global_high_freq_vars - jpn_high_freq_vars])
    both_in_jpn_and_global = np.array([(jpn_vars.get(k, 0), global_vars.get(k, 0)) for k in global_high_freq_vars & jpn_high_freq_vars])

    ax.axvspan(0.5, 1.05, color='lightgray', alpha=0.6, lw=0)
    ax.axhspan(0.5, 1.05, color='lightgray', alpha=0.6, lw=0)

    ax.scatter(both_in_jpn_and_global.T[0], both_in_jpn_and_global.T[1], marker='o', ec='darkgray', fc='white', label="AF > 0.5 \nin both 54KJPN and 1KGP")
    ax.scatter(only_in_jpn.T[0], only_in_jpn.T[1], marker='v', color='tab:orange', label="AF > 0.5 \nonly in 54KJPN")
    ax.scatter(only_in_global.T[0], only_in_global.T[1], marker='^', color='tab:gray', label="AF > 0.5 \nonly in 1KGP")

    ax.legend(bbox_to_anchor=(0.25, 0.25), loc='center', fontsize=11)

    ax.set_xlim(-0.05, 1.05)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels([0, 0.25, 0.5, 0.75, 1], size=16)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticklabels([0, 0.25, 0.5, 0.75, 1], size=16)
    ax.set_ylim(-0.05, 1.05)
    ax.set_aspect('equal')
    ax.set_xlabel("Allele Freq. in 54KJPN", size=20)
    ax.set_ylabel("Allele Freq. in 1KGP", size=20)
    fig.tight_layout()
    fig.savefig("./figures/4b_jpn_vs_global.pdf")

    notable_vars = [k for k in jpn_high_freq_vars - global_high_freq_vars if global_vars.get(k, 0) < 0.05]
    notable_vars.sort(key=lambda v: jpn_vars[v], reverse=True)

    for v in [(k, jpn_vars[k], global_vars.get(k, 0)) for k in notable_vars]:
        text = v[0].replace("_human", "").upper()
        print(text, v[1], v[2])

if __name__ == '__main__':
    analyze_jpn_vs_global()