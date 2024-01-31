#!/usr/bin/env python3
from misc import AA2COLOR
import gpcrdb
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = "Arial"
from typing import List, Dict
import os
from utils import VariationType, Segment
import json
import numpy as np
import sys


def plot_G_protein_common_contact_positions(vars: Dict):
    fig, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=300, sharex=True, sharey=True)
    af_values = np.array([[vars[receptor][gnum][0] for gnum in vars[receptor].keys()] for receptor in vars.keys()]).flatten()
    bins = np.insert(np.linspace(1 / 50, 1, 50), 0, [0, sys.float_info.epsilon])
    ax.hist(af_values, bins=bins, histtype='step', lw=1, color='tab:gray')
    ax.hist(af_values, bins=bins, color='tab:gray')
    for receptor in vars.keys():
        for gnum in vars[receptor].keys():
            af = vars[receptor][gnum][0]
            anno = vars[receptor][gnum][1]
            if af > 0.05:
                x, y = af, 1
                ha = 'left'
                if af < 0.2:
                    dx, dy = 0.05, 10
                    ax.arrow(af, 1, dx=dx, dy=dy, lw=0.5)
                    x += dx
                    y += dy
                    ha = 'center'
                elif 0.2 < af < 0.8:
                    ha = 'center'
                else:
                    ha = 'right'
                ax.text(x, y, anno, va='bottom', ha=ha)
    ax.set_xlabel("Allele Freq.")
    ax.set_yscale('log')
    ax.set_ylabel("Number of Variants")
    fig.tight_layout()
    fig.savefig("G_protein_common_contact_positions.pdf")

def analyze_G_protein_common_contact_positions() -> Dict:
    # See Fig. 3c for detail
    # https://doi.org/10.1038/s41467-022-34055-5
    common_contacts = {"3.50", "3.53", "3.54", "34.50", "34.51", "34.55", "5.65", "5.68", "6.32", "6.33", "6.36", "6.37", "7.56", "8.47"}

    vars = {}
    for receptor in gpcrdb.get_filtered_receptor_list("receptors.json"):
        entry_name = receptor['entry_name']
        receptor_class = receptor['receptor_class']
        if receptor_class != 'Class A (Rhodopsin)':
            continue

        dpath = os.path.join(receptor_class, entry_name)
        with open(os.path.join(dpath, '38KJPN-CDS.csv')) as f:
            gnum_af = {gnum: [0, ""] for gnum in common_contacts}
            for l in f.readlines():
                cols = l.strip().split('\t')
                if l.startswith('#'):
                    header = cols
                    continue
                af = float(cols[header.index('AF')])
                gnum = cols[header.index('Generic_Num')]
                vtype = VariationType.name_of(cols[header.index('#Var_Type')])
                if vtype == VariationType.MISSENSE and gnum in common_contacts:
                    ref_aa, res_num, alt_aa = cols[header.index('Ref_AA')], cols[header.index('Res_Num')], cols[header.index('Alt_AA')]
                    annotation = "{}\n{}{}$^{{{}}}${}\nAF = {:.2f}".format(entry_name.replace("_human", "").upper(), ref_aa, res_num, gnum, alt_aa, af)
                    gnum_af[gnum] = [af, annotation]
            vars[entry_name] = gnum_af
    return vars

def plot_G_protein_contact_positions(af: Dict, anno: List):
    titles = ['Gs complex', 'Gi complex', 'Gq complex', 'common']
    ylabels = ['Gs', 'Gi/o', 'Gq/11', 'promiscuous']

    bins = np.insert(np.linspace(1 / 20, 1, 20), 0, [0, sys.float_info.epsilon])

    fig, axes = plt.subplots(len(titles), len(ylabels), figsize=(8, 8), dpi=300, sharex=True)
    for y, primary_coupled_G_protein in enumerate(ylabels):
        for x, positions in enumerate(titles):
            ax = axes[y][x]
            if x == 0:
                ax.set_ylabel(primary_coupled_G_protein, size=12)
            if y == 0:
                ax.set_title(positions, size=12)
            if y == 3:
                ax.set_xlabel("Allele Freq.", size=10)
            
            ax.hist(af[primary_coupled_G_protein][positions], bins=bins, color='tab:gray')

            for a in anno:
                if a[0] == primary_coupled_G_protein and a[1] == positions:
                    ax.text(a[2], 1, a[3], size=10, va='bottom', ha='center' if a[2] > 0.5 else 'left')

    fig.text(0.5, 0.95, "G protein contact positions", size=14, ha='center', va='top')
    fig.text(0.05, 0.5, "Number of missense receptor variants by primary G protein", rotation=90, size=14, ha='right', va='center')
    
    #fig.tight_layout()
    fig.savefig("G_protein_contact_positions.pdf")

def analyze_G_protein_contact_positions() -> Dict:
    # See Fig. 3c and Sup. Fig. 2a-c for detail
    # https://doi.org/10.1038/s41467-022-34055-5
    common_contacts = {"3.50", "3.53", "3.54", "34.50", "34.51", "34.55", "5.65", "5.68", "6.32", "6.33", "6.36", "6.37", "7.56", "8.47"}
    gs_contacts = {"3.54", "3.55", "34.51", "34.54", "34.55", "5.64", "5.68", "5.69", "5.71", "5.72", "5.74", "5.75", "5.77", "8.48"}
    gi_contacts = {"12.49", "2.40", "3.50", "3.53", "34.52", "34.55", "5.71", "6.32", "7.56", "8.47", "8.49"}
    gq_contacts = {"2.37", "2.39", "2.40", "3.49", "34.51", "34.53", "34.55", "34.56", "34.57", "4.38", "4.39", "6.30", "6.33", "8.48", "8.49"}
    gs_specific_contacts = gs_contacts - common_contacts
    gi_specific_contacts = gi_contacts - common_contacts
    gq_specific_contacts = gq_contacts - common_contacts
    contact_positions = {
        "Gs complex": gs_specific_contacts, 
        "Gi complex": gi_specific_contacts, 
        "Gq complex": gq_specific_contacts, 
        "common": common_contacts
    }
    print(contact_positions)

    uniprot2dpath = {receptor['entry_name'].replace('_human', '').upper(): os.path.join(receptor['receptor_class'], receptor['entry_name']) 
                     for receptor in gpcrdb.get_filtered_receptor_list("receptors.json")}

    ret = {}
    annotations = []
    couplings = gpcrdb.primary_coupled_receptors()
    for primary_coupled_G_protein in couplings.keys():
        ret[primary_coupled_G_protein] = {}
        for uniprot_id in couplings[primary_coupled_G_protein]:
            dpath = uniprot2dpath[uniprot_id]
            for contact_type in contact_positions.keys():
                positions = contact_positions[contact_type]
                with open(os.path.join(dpath, '38KJPN-CDS.csv')) as f:
                    for l in f.readlines():
                        cols = l.strip().split('\t')
                        if l.startswith('#'):
                            header = cols
                            continue
                        af = float(cols[header.index('AF')])
                        gnum = cols[header.index('Generic_Num')]
                        vtype = VariationType.name_of(cols[header.index('#Var_Type')])
                        if vtype == VariationType.MISSENSE and gnum in positions:
                            ref_aa, res_num, alt_aa = cols[header.index('Ref_AA')], cols[header.index('Res_Num')], cols[header.index('Alt_AA')]
                            annotation = "{}\n{}{}$^{{{}}}${}\nAF = {:.2f}".format(uniprot_id, ref_aa, res_num, gnum, alt_aa, af)
                            if af > 0.2:
                                annotations.append([primary_coupled_G_protein, contact_type, af, annotation])
                            arr = ret[primary_coupled_G_protein].get(contact_type, [])
                            arr.append(af)
                            ret[primary_coupled_G_protein][contact_type] = arr
    print(annotations)
    return ret, annotations

def main():
    # plot_gene_map()
    # calls = analyze_calls()
    # plot_calls(*calls)
    # var_type = analyze_var_type()
    # plot_var_type(*var_type)
    # var_seg = analyze_var_seg()
    # plot_var_seg(*var_seg)
    # plot_var_seg_percent(*var_seg)
    # af = analyze_allele_freq()
    # plot_allele_freq(af)
    # nter = analyze_Nter()
    # plot_Nter(nter)
    # cter = analyze_Cter()
    # plot_cter(cter)
    # vars = analyze_high_freq_vars()
    # plot_high_freq_vars(vars)
    # pos = analyze_family_A_pos()
    # plot_family_A_pos(*pos)
    # aa_stats, codon_stats = analyze_arginine_3x50()
    # plot_arginine_3x50(aa_stats, codon_stats)
    vars = analyze_G_protein_common_contact_positions()
    plot_G_protein_common_contact_positions(vars)
    # af, anno = analyze_G_protein_contact_positions()
    # plot_G_protein_contact_positions(af, anno)
    pass

if __name__ == '__main__':
    main()