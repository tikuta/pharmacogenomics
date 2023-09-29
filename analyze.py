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




def plot_allele_freq(af: Dict):
    fig, ax = plt.subplots(1, 1, figsize=(6, 2), dpi=300)

    ax.axhspan(-0.5, len(Segment) - 0.5, xmin=0.5, color='whitesmoke')
    ax.set_ylim(-0.5, len(Segment) - 0.5)
    for i, seg in enumerate(Segment):
        ax.scatter(af[seg.value], [i] * len(af[seg.value]), color=seg.color)

    ax.set_yticks(range(len(Segment)))
    ax.set_yticklabels([seg.value for seg in Segment])
    ax.invert_yaxis()
    ax.set_xlabel("Allele Freq.")
    fig.tight_layout()
    fig.savefig("allele_freq.pdf")

def analyze_allele_freq() -> Dict:
    af = {}
    for receptor in gpcrdb.get_filtered_receptor_list("receptors.json"):
        entry_name = receptor['entry_name']
        receptor_class = receptor['receptor_class']
        dpath = os.path.join(receptor_class, entry_name)
        with open(os.path.join(dpath, '38KJPN-CDS.csv')) as f:
            for l in f.readlines():
                if l.startswith('#'):
                    continue
                cols = l.split('\t')
                var_type = cols[0]
                seg = cols[10]
                freq = float(cols[14])
                if seg == 'None':
                    print(dpath, l)
                if var_type == VariationType.MISSENSE.name:
                    freqs = af.get(seg, [])
                    freqs.append(freq)
                    af[seg] = freqs
    return af

def plot_cter(cter: Dict):
    fig, ax = plt.subplots(1, 1, figsize=(6, 2), dpi=300)

    ax.set_facecolor('whitesmoke')
    for i, k in enumerate(cter.keys()):
        if k == "non-S/T → non-S/T":
            continue
        freqs = cter[k]['freqs']
        labels = cter[k]['labels']
        ax.scatter(freqs, [i] * len(freqs), color=Segment.Cterm.color)
        
        pos = []
        for f, l in zip(freqs, labels):
            if f > 0.5:
                x = f
                ha = 'center'
                if i == 1:
                    y = i + 0.1
                    va = 'bottom'
                    for p in pos:
                        if y != p[1]:
                            continue
                        if abs(x - p[0]) < 0.1:
                            va = 'top'
                            y -= 0.3
                            break
                else:
                    y, va = i - 0.2, 'top'
                    if l.startswith('NMUR2'):
                        x += 0.005
                        ha = 'left'
                    elif l.startswith('RXFP4'):
                        x -= 0.005
                        ha = 'right'
                ax.text(x, y, l, va=va, ha=ha)
                pos.append([x, y])
                print(l, f)

    ax.set_yticks(range(len(cter)))
    ax.set_yticklabels(list(cter.keys()))
    ax.set_ylim(-1, 2)
    ax.set_ylabel("Ser/Thr")
    ax.set_xlabel("Allele Freq.")
    ax.set_xlim(left=0.5, right=1.09)
    fig.tight_layout()
    fig.savefig("cter.pdf")

def analyze_Cter() -> Dict:
    cter = {}
    for receptor in gpcrdb.get_filtered_receptor_list("receptors.json"):
        entry_name = receptor['entry_name']
        receptor_class = receptor['receptor_class']
        dpath = os.path.join(receptor_class, entry_name)
        with open(os.path.join(dpath, 'ensembl.json')) as f:
            display_name = json.load(f)['display_name']
        with open(os.path.join(dpath, '38KJPN-CDS.csv')) as f:
            for l in f.readlines():
                if l.startswith('#'):
                    continue
                cols = l.split('\t')
                var_type = cols[0]
                seg = cols[10]
                if seg != 'C-term' or var_type != VariationType.MISSENSE.name:
                    continue
                aa_ref, aa_alt = cols[6], cols[9]
                freq = float(cols[14])
                res_num = cols[3]
                ser_and_thr = ('S', 'T')
                if aa_ref in ser_and_thr:
                    if aa_alt in ser_and_thr:
                        key = "S/T → S/T"
                        continue
                    else:
                        key = "Loss"
                else:
                    if aa_alt in ser_and_thr:
                        key = "Gain"
                    else:
                        key = "non-S/T → non-S/T"
                        continue
                d = cter.get(key, {"freqs": [], "labels": []})
                d["freqs"].append(freq)
                d["labels"].append(f"{display_name}\n{aa_ref}{res_num}{aa_alt}")
                cter[key] = d
    return cter

def plot_Nter(glyco: Dict):
    fig, ax = plt.subplots(1, 1, figsize=(6, 2.5), dpi=300)

    ax.set_facecolor('whitesmoke')
    for i, k in enumerate(glyco.keys()):
        n_labels = [v[0] for v in glyco[k] if v[2]]
        n_freqs = [v[1] for v in glyco[k] if v[2]]
        ax.scatter(n_freqs, [i] * len(n_freqs), facecolor='white', edgecolors=Segment.Nterm.color, zorder=10, label="N-glycosylation motif (N-X-S/T)" if i == 0 else None)
        for f, l in zip(n_freqs, n_labels):
            if f > 0.5:
                ax.text(f, i - 0.2 if i == 0 else i + 0.2, l, va='top' if i == 0 else 'bottom', ha='center')
        
        o_labels = [v[0] for v in glyco[k] if not v[2]]
        o_freqs = [v[1] for v in glyco[k] if not v[2]]
        ax.scatter(o_freqs, [i] * len(o_freqs), facecolor=Segment.Nterm.color, zorder=1, label="O-glycosylation (S/T)" if i == 0 else None)
    
    ax.legend(bbox_to_anchor=(0.5, 1), loc='lower center', ncol=2)
        
    #ax.legend(bbox_to_anchor=(0.5, 1.1), loc='lower center', ncol=2)
    ax.set_xlabel("Allele Freq.")
    ax.set_xlim(left=0.5)
    ax.set_ylim(-1, 2)
    ax.set_yticks(range(len(glyco)))
    ax.set_yticklabels(glyco.keys())
    ax.set_ylabel("Asn/Ser/Thr")
    fig.tight_layout()
    
    fig.savefig("nter.pdf")


def analyze_Nter() -> Dict:
    glyco = {}
    for receptor in gpcrdb.get_filtered_receptor_list("receptors.json"):
        entry_name = receptor['entry_name']
        receptor_class = receptor['receptor_class']
        dpath = os.path.join(receptor_class, entry_name)
        with open(os.path.join(dpath, 'ensembl.json')) as f:
            display_name = json.load(f)['display_name']
        with open(os.path.join(dpath, '38KJPN-CDS.csv')) as f:
            with open(os.path.join(dpath, "match.json")) as j:
                seq = ''.join([r['ensembl_amino_acid'] for r in json.load(j)['residues']])
            for l in f.readlines():
                if l.startswith('#'):
                    continue
                cols = l.split('\t')
                var_type = cols[0]
                res_num = int(cols[3])
                aa_ref = cols[6]
                aa_alt = cols[9]
                seg = cols[10]
                if seg != Segment.Nterm.value or var_type != VariationType.MISSENSE.name:
                    continue
                freq = float(cols[14])
                SER_THR = ('S', 'T')
                ASN = ('N', )
                
                label = f"{display_name}\n{aa_ref}{res_num}{aa_alt}"
                
                if aa_ref in SER_THR and aa_alt not in SER_THR:
                    key = "Loss"
                elif aa_ref not in SER_THR and aa_alt in SER_THR:
                    key = "Gain"

                is_related_to_n_glyco = False
                if aa_ref in ASN or aa_alt in ASN:
                    assert(seq[res_num - 1] == aa_ref)
                    if seq[res_num] != 'P' and seq[res_num + 1] in SER_THR: # [N]-X-S/T
                        is_related_to_n_glyco = True
                        #label += " " + seq[res_num - 1: res_num + 2]
                        if aa_ref in ASN:
                            key = "Loss"
                        else:
                            key = "Gain"
                    else: # not related with [N]-X-S/T
                        if key: # N <-> S/T (simple S/T gain/loss)
                            pass
                        else: # N <-> non-S/T
                            continue
                else:
                    if key: # non-N <-> non-S/T
                        pass
                    else: # non-N <-> non-S/T
                        continue

                a = glyco.get(key, [])
                a.append([label, freq, is_related_to_n_glyco])
                glyco[key] = a
    return glyco

def plot_high_freq_vars(vars: List):
    fig, ax = plt.subplots(1, 1, figsize=(4, 12), dpi=300)
    ax.set_facecolor('whitesmoke')

    vars.sort(key=lambda v: (v[1], v[3]))

    for i, var in enumerate(vars):
        ax.barh(i, var[1], color=var[2])
    
    A_ticks, A_ticklabels = [], []
    non_A_ticks, non_A_ticklabels = [], []
    for i, v in enumerate(vars):
        if v[3]:
            A_ticks.append(i)
            A_ticklabels.append(v[0])
        else:
            non_A_ticks.append(i)
            non_A_ticklabels.append(v[0])
    ax.set_yticks(A_ticks)
    ax.set_yticklabels(A_ticklabels, size=7)
    ax.set_yticks(non_A_ticks, minor=True)
    ax.set_yticklabels(non_A_ticklabels, size=7, minor=True, color='tab:gray')
    ax.set_xlabel("Allele Freq.")
    ax.set_ylim(-0.8, len(vars) - 0.3)
    ax.set_xticks([0.5, 0.6, 0.8, 1.0])
    ax.set_xticklabels([0.5, 0.6, 0.8, 1.0])
    ax.set_xlim(left=0.5)
    fig.tight_layout()
    fig.savefig("high_freq_vars.pdf")


def analyze_high_freq_vars() -> List:
    vars = []
    labels = ["N-term", "TM1", "ICL1", "TM2", "ECL1", "TM3", "ICL2", "TM4", "ECL2", "TM5", "ICL3", "TM6", "ECL3", "TM7", "H8", "ICL4", "C-term", "None"]
    cmap = plt.get_cmap('rainbow', len(labels) - 1)
    colors = [cmap(i) for i in range(len(labels) - 1)] + ['tab:gray']
    for receptor in gpcrdb.get_filtered_receptor_list("receptors.json"):
        entry_name = receptor['entry_name']
        receptor_class = receptor['receptor_class']
        dpath = os.path.join(receptor_class, entry_name)
        with open(os.path.join(dpath, 'ensembl.json')) as f:
            display_name = json.load(f)['display_name']
        with open(os.path.join(dpath, '38KJPN-CDS.csv')) as f:
            for l in f.readlines():
                if l.startswith('#'):
                    continue
                cols = l.split('\t')
                var_type = cols[0]
                res_num = cols[3]
                aa_ref = cols[6]
                aa_alt = cols[9]
                seg = cols[10]
                gen_num = cols[11]
                freq = float(cols[14])
                if seg == 'None':
                    print(dpath, l)
                if var_type == VariationType.MISSENSE.name:
                    if freq > 0.5:
                        label = "{} {}{}{} ({})".format(display_name, aa_ref, res_num, aa_alt, gen_num if gen_num != 'None' else seg)
                        color = colors[labels.index(seg)]
                        classA = True if receptor_class == 'Class A (Rhodopsin)' else False
                        var = [label, freq, color, classA]
                        vars.append(var)
    print(sorted([v[0] for v in vars if v[3]], key=lambda v: v.split('(')[1]))
    return vars

def plot_family_A_pos(count: Dict, g_num_count: Dict):
    g_nums = list(g_num_count.keys())
    g_nums.sort()

    fig, ax = plt.subplots(1, 1, figsize=(3, 6), dpi=300)
    ax.invert_yaxis()
    for i, g_num in enumerate(g_nums):
        color = Segment.generic_number_of(g_num).color
        ax.barh(i, g_num_count[g_num], color='lightgray', height=1, zorder=-200)
        ax.barh(i, count.get(g_num, 0), color=color, height=1, zorder=-100)

    ax.set_yticks([g_nums.index(g_num) for g_num in g_nums if g_num.endswith('.50')])
    ax.set_yticklabels([g_num for g_num in g_nums if g_num.endswith('.50')], size=8)
    ax.set_ylabel("Generic number (BW number)")
    ax.set_xlim(0, max(g_num_count.values()))
    ax.set_xlabel("Number of family A GPCRs")
    ax.set_ylim(len(g_nums), -1)
    ax.set_xticks([0, 50, 100, 150, 200, 250, max(g_num_count.values())])
    ax.set_xticklabels([0, 50, 100, 150, 200, 250, max(g_num_count.values())])
    # ax.legend(ncol=2, bbox_to_anchor=(0.5, 1.1), loc='lower center')
    fig.tight_layout()
    fig.savefig("family_A_pos.pdf")


def analyze_family_A_pos() -> List[Dict]:
    count = {}
    g_num_count = {}
    for receptor in gpcrdb.get_filtered_receptor_list("receptors.json"):
        entry_name = receptor['entry_name']
        receptor_class = receptor['receptor_class']
        if receptor_class != 'Class A (Rhodopsin)':
            continue

        dpath = os.path.join(receptor_class, entry_name)
        missense, silent, nonsense = set(), set(), set()
        with open(os.path.join(dpath, '38KJPN-CDS.csv')) as f:
            for l in f.readlines():
                if l.startswith('#'):
                    continue
                cols = l.split('\t')
                var_type = VariationType.name_of(cols[0])
                # res_num = cols[3]
                # aa_ref = cols[6]
                # aa_alt = cols[9]
                # seg = cols[10]
                gen_num = cols[11]
                # freq = float(cols[14])
                if gen_num != 'None':
                    if var_type == VariationType.MISSENSE:
                        missense.add(gen_num)
                    elif var_type == VariationType.SILENT:
                        silent.add(gen_num)
                    elif var_type == VariationType.NONSENSE:
                        nonsense.add(gen_num)
        for g_num in missense:
            count[g_num] = count.get(g_num, 0) + 1

        with open(os.path.join(dpath, 'match.json')) as f:
            for r in json.load(f)['residues']:
                gen_num = r['generic_number']
                if gen_num:
                    g_num_count[gen_num] = g_num_count.get(gen_num, 0) + 1

    return count, g_num_count

def plot_arginine_3x50(aa_stats: Dict, codon_stats: Dict):
    fig, axes = plt.subplots(2, 1, figsize=(8, 2), dpi=300, sharex=True)

    ax = axes[0]
    aas = list(aa_stats.keys())
    aas.sort(key=lambda aa: (aa != 'Unassigned', aa_stats[aa]))
    left = 0
    for aa in aas[::-1]:
        delta = aa_stats[aa]
        print(aa, delta)
        ax.barh(0, delta, left=left, color=AA2COLOR.get(aa, 'tab:gray'), linewidth=0.5, edgecolor='k')
        if aa == "R":
            ax.text(left + delta / 2, 0, "{}\n({})".format(aa, delta), ha='center', va='center')
        left += delta
    ax.set_axis_off()

    ax = axes[1]
    codons = list(codon_stats.keys())
    codons.sort(key=lambda codon: codon_stats[codon])
    left = 0
    for codon in codons[::-1]:
        delta = codon_stats[codon]
        print(codon, delta)
        ax.barh(0, delta, left=left, linewidth=0.5, edgecolor='k')
        y, va = 0.45 if delta < 20 else 0, 'bottom' if delta < 20 else 'center'
        ax.text(left + delta / 2, 0, "{}\n({})".format(codon, delta), ha='center', va='center')
        left += delta
    ax.set_axis_off()

    fig.tight_layout()
    fig.savefig("arginine_3x50.pdf")

def analyze_arginine_3x50() -> List[Dict]:
    aa_stats = {}
    codon_stats = {}
    for receptor in gpcrdb.get_filtered_receptor_list("receptors.json"):
        entry_name = receptor['entry_name']
        receptor_class = receptor['receptor_class']
        if receptor_class != 'Class A (Rhodopsin)':
            continue

        dpath = os.path.join(receptor_class, entry_name)
        with open(os.path.join(dpath, 'match.json')) as f:
            residues = json.load(f)['residues']
            g_nums = {r['generic_number']: r['ensembl_amino_acid'] for r in residues if r['generic_number'] is not None}
            res_nums = {r['generic_number']: r['ensembl_sequence_number'] for r in residues if r['generic_number'] is not None}

            roi = "3.50"
            aa = g_nums.get(roi, "Unassigned")
            aa_stats[aa] = aa_stats.get(aa, 0) + 1

            if aa == 'R':
                with open(os.path.join(dpath, 'cds.json')) as f:
                    cds = json.load(f)
                res_num = res_nums[roi]
                codon = cds['nucleotides'][res_num * 3 - 3: res_num * 3]
                codon_stats[codon] = codon_stats.get(codon, 0) + 1
    return aa_stats, codon_stats

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