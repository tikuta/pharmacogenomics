#!/usr/bin/env python3
from misc import AA2COLOR
import gpcrdb
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = "Arial"
from typing import List, Dict
import os
from utils import VariationType
import json
from misc import Segment
import numpy as np

def plot_calls(num_cds: int, num_gene: int):
    fig, ax = plt.subplots(1, 1, figsize=(4, 2), dpi=300)

    ax.barh(0, num_cds, color='tab:orange')
    coding_text = "Coding region\n{:,} calls ({:.1f}%)".format(num_cds, num_cds / num_gene * 100)
    ax.text(num_cds / 2, -0.45, coding_text, ha='center', va='top', size=12)

    ax.barh(0, num_gene - num_cds, left=num_cds, color='tab:gray', alpha=0.6)
    non_coding_text = "Non-coding region\n{:,} calls ({:.1f}%)".format(num_gene - num_cds, (num_gene - num_cds) / num_gene * 100)
    ax.text((num_gene - num_cds) / 2, 0, non_coding_text, ha='center', va='center', size=12)

    ax.set_xlim(0, num_gene)
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig("calls.pdf")

def analyze_calls() -> List[int]:
    num_cds, num_gene = 0, 0
    for receptor in gpcrdb.get_filtered_receptor_list("receptors.json"):
        entry_name = receptor['entry_name']
        receptor_class = receptor['receptor_class']
        dpath = os.path.join(receptor_class, entry_name)

        calls_gene = set()
        with open(os.path.join(dpath, '38KJPN-gene.vcf')) as f:
            for l in f.readlines():
                cols = l.split('\t')
                chromosome = cols[0]
                position = int(cols[1])
                calls_gene.add((chromosome, position))

        calls_cds = set()
        with open(os.path.join(dpath, '38KJPN-CDS.vcf')) as f:
            for l in f.readlines():
                cols = l.split('\t')
                chromosome = cols[0]
                position = int(cols[1])
                calls_cds.add((chromosome, position))

        assert(calls_cds.issubset(calls_gene))
        num_cds += len(calls_cds)
        num_gene += len(calls_gene)
    print("Total", num_cds, num_gene)
    return num_cds, num_gene

def plot_var_type(num_missense: int, num_silent: int, num_nonsense: int):
    total = num_missense + num_silent + num_nonsense

    fig, ax = plt.subplots(1, 1, figsize=(4, 2), dpi=300)

    ax.barh(0, num_missense, color='tab:orange')
    missense_text = "Missense\n{:,} SNVs\n({:.1f}%)".format(num_missense, num_missense / total * 100)
    ax.text(num_missense / 2, 0, missense_text, ha='center', va='center', size=12)

    ax.barh(0, num_silent, left=num_missense, color='tab:gray', alpha=0.6)
    silent_text = "Silent\n{:,} SNVs\n({:.1f}%)".format(num_silent, num_silent / total * 100)
    ax.text(num_missense + num_silent / 2, 0, silent_text, ha='center', va='center', size=12)

    ax.barh(0, num_nonsense, left=num_missense + num_silent, color='tab:gray')
    nonsense_text = "Nonsense\n{:,} SNVs\n({:.1f}%)".format(num_nonsense, num_nonsense / total * 100)
    ax.text(num_missense + num_silent + num_nonsense / 2, -0.45, nonsense_text, ha='center', va='top', size=12)

    ax.set_xlim(0, total)
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig("var_type.pdf")

def analyze_var_type():
    num_missense, num_silent, num_nonsense = 0, 0, 0
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
                if var_type == VariationType.MISSENSE.name:
                    num_missense += 1
                    continue
                elif var_type == VariationType.SILENT.name:
                    num_silent += 1
                    continue
                elif var_type == VariationType.NONSENSE.name:
                    num_nonsense += 1
                    continue
                raise Exception("Illegal variation type `{}`".format(var_type))
    print("Total", num_missense, num_silent, num_nonsense)
    return num_missense, num_silent, num_nonsense

def plot_var_seg(seg_missense: Dict, seg_silent: Dict, seg_nonsense: Dict):
    # Ignore ICL4 or None (segment unassigned)
    xlabels = ["N-term", "TM1", "ICL1", "TM2", "ECL1", "TM3", "ICL2", "TM4", "ECL2", "TM5", "ICL3", "TM6", "ECL3", "TM7", "H8", "C-term"]
    fig, ax = plt.subplots(1, 1, figsize=(4, 4), dpi=300)

    for x, xlabel in enumerate(xlabels):
        ax.bar(x - 0.3, seg_missense[xlabel], width=0.3, color='tab:orange')
        ax.bar(x, seg_silent[xlabel], width=0.3, color='tab:gray', alpha=0.6)
        ax.bar(x + 0.3, seg_nonsense[xlabel], width=0.3, color='tab:gray')

    ax.set_xticks(range(len(xlabels)))
    ax.set_xticklabels(xlabels, rotation=90)
    ax.set_ylabel("Number of SNVs")
    fig.tight_layout()
    fig.savefig("var_seg.pdf")

def plot_var_seg_percent(seg_missense: Dict, seg_silent: Dict, seg_nonsense: Dict):
    xlabels = ["N-term", "TM1", "ICL1", "TM2", "ECL1", "TM3", "ICL2", "TM4", "ECL2", "TM5", "ICL3", "TM6", "ECL3", "TM7", "H8", "ICL4", "C-term", "None"]
    cmap = plt.get_cmap('rainbow', len(xlabels) - 1)
    colors = [cmap(i) for i in range(len(xlabels) - 1)] + ['tab:gray']

    total_missense = sum([v for v in seg_missense.values()])
    total_silent = sum([v for v in seg_silent.values()])
    total_nonsense = sum([v for v in seg_nonsense.values()])
    total = total_missense + total_silent + total_nonsense

    fig, ax = plt.subplots(1, 1, figsize=(4, 2), dpi=300)

    left = 0
    for x in range(len(xlabels)):
        width = seg_missense.get(xlabels[x], 0) / total_missense * 100
        ax.barh(2, width, height=0.7, left=left, color=colors[x], edgecolor='k', lw=0.2)
        
        label = xlabels[x]
        if label.startswith('TM') or label.endswith('-term'):
            ax.text(left + width / 2, 2.4, label, ha='center', va='bottom', size=6)
        left += width

    left = 0
    for x in range(len(xlabels)):
        width = seg_silent.get(xlabels[x], 0) / total_silent * 100
        ax.barh(1, width, height=0.7, left=left, color=colors[x], edgecolor='k', lw=0.2)
        left += width

    left = 0
    for x in range(len(xlabels)):
        width = seg_nonsense.get(xlabels[x], 0) / total_nonsense * 100
        ax.barh(0, width, height=0.7, left=left, color=colors[x], edgecolor='k', lw=0.2)
        left += width

    ax.set_ylim(-0.5, 2.7)
    ax.set_yticks(range(3))
    ax.set_yticklabels(["Nonsense", "Silent", "Missense"])
    ax.set_xlabel("Fraction [%]")
    fig.tight_layout()
    fig.savefig("var_seg_percent.pdf")

def analyze_var_seg():
    seg_missense, seg_silent, seg_nonsense = {}, {}, {}
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
                if seg == 'None':
                    print(dpath, l)
                if var_type == VariationType.MISSENSE.name:
                    seg_missense[seg] = seg_missense.get(seg, 0) + 1
                elif var_type == VariationType.SILENT.name:
                    seg_silent[seg] = seg_silent.get(seg, 0) + 1
                elif var_type == VariationType.NONSENSE.name:
                    seg_nonsense[seg] = seg_nonsense.get(seg, 0) + 1
    print(seg_missense)
    print(seg_silent)
    print(seg_nonsense)
    return seg_missense, seg_silent, seg_nonsense

def plot_allele_freq(af: Dict):
    fig, ax = plt.subplots(1, 1, figsize=(6, 4), dpi=300)

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

    fig, ax = plt.subplots(1, 1, figsize=(6, 4), dpi=300)
    ax2 = ax.twinx()
    min_receptor = 100
    ax.axhspan(0, min_receptor, color='whitesmoke')
    for i, g_num in enumerate(g_nums):
        color = Segment.generic_number_of(g_num).color
        ax.bar(i, g_num_count[g_num], color='lightgray')
        ax.bar(i, count.get(g_num, 0), color=color)
        if g_num_count[g_num] >= min_receptor: # at least 10% of receptors has the generic number
            percentage = count.get(g_num, 0) / g_num_count[g_num] * 100
            ax2.scatter(i, percentage, color=color, marker='s', s=0.2)
            if percentage > 20:
                ax2.text(i, percentage, "{}\n{} GPCRs ({:.1f}%)".format(g_num, count[g_num], percentage), va='bottom', ha='center', size=8)

    ax.set_xticks([g_nums.index(g_num) for g_num in g_nums if g_num.endswith('.50')])
    ax.set_xticklabels([g_num for g_num in g_nums if g_num.endswith('.50')], rotation=90, size=8)
    ax.set_xlabel("Generic number (BW number)")
    ax.set_ylim(0, max(g_num_count.values()))
    ax.set_ylabel("Number of family A GPCRs")
    ax.set_yticks([0, 50, 100, 150, 200, 250, max(g_num_count.values())])
    ax.set_yticklabels([0, 50, 100, 150, 200, 250, max(g_num_count.values())])
    ax.set_yticks([min_receptor], minor=True)
    ax.set_yticklabels([min_receptor], minor=True)
    ax2.set_ylabel("GPCRs that can carry missense variations [%]")
    ax2.set_ylim(0, 35)
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

def plot_amino_acid(amino_acids: Dict):
    fig, ax = plt.subplots(1, 1, figsize=(6, 2), dpi=300)
    keys = list(amino_acids.keys())[::-1]
    for i, pos in enumerate(keys):
        aas = list(amino_acids[pos].keys())
        aas.sort(key=lambda aa: (aa != 'Unassigned', amino_acids[pos][aa]))
        left = 0
        for aa in aas[::-1]:
            delta = amino_acids[pos][aa]
            print(pos, aa, delta)
            ax.barh(i, delta, left=left, color=AA2COLOR.get(aa, 'tab:gray'), linewidth=0.5, edgecolor='k')
            ax.text(left + delta / 2, i, aa, ha='center', va='center')
            left += delta
    ax.set_xticks([0, 50, 100, 150, 200, 250, left])
    ax.set_xticklabels([0, 50, 100, 150, 200, 250, left])
    ax.set_xlim(0, left)
    ax.set_xlabel("Number of GPCRs")
    ax.set_yticks(range(len(keys)))
    ax.set_yticklabels(keys)
    ax.set_ylabel("Generic number")
    fig.tight_layout()
    fig.savefig("amino_acid.pdf")

def analyze_amino_acid() -> Dict:
    amino_acids = {'3.50': {}, '8.51': {}, '34.57': {}}
    for receptor in gpcrdb.get_filtered_receptor_list("receptors.json"):
        entry_name = receptor['entry_name']
        receptor_class = receptor['receptor_class']
        if receptor_class != 'Class A (Rhodopsin)':
            continue

        dpath = os.path.join(receptor_class, entry_name)
        with open(os.path.join(dpath, 'match.json')) as f:
            residues = json.load(f)['residues']
            g_nums = {r['generic_number']: r['ensembl_amino_acid'] for r in residues if r['generic_number'] is not None}
            for g_num in amino_acids.keys():
                d = amino_acids[g_num]
                if g_num in g_nums.keys():
                    aa = g_nums[g_num]
                else:
                    aa = "Unassigned"
                d[aa] = d.get(aa, 0) + 1
                amino_acids[g_num] = d
    return amino_acids

def main():
    # plot_gene_map()
    # calls = analyze_calls()
    # plot_calls(*calls)
    # var_type = analyze_var_type()
    # plot_var_type(*var_type)
    # var_seg = analyze_var_seg()
    # plot_var_seg(*var_seg)
    # plot_var_seg_percent(*var_seg)
    af = analyze_allele_freq()
    plot_allele_freq(af)
    # nter = analyze_Nter()
    # plot_Nter(nter)
    # cter = analyze_Cter()
    # plot_cter(cter)
    # vars = analyze_high_freq_vars()
    # plot_high_freq_vars(vars)
    # pos = analyze_family_A_pos()
    # plot_family_A_pos(*pos)
    # aa = analyze_amino_acid()
    # plot_amino_acid(aa)
    pass

if __name__ == '__main__':
    main()