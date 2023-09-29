#!/usr/bin/env python3
from misc import VCF_JPN_GENE_FILENAME, VCF_JPN_CDS_FILENAME, CSV_JPN_CDS_FILENAME
import gpcrdb
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = "Arial"
import os
from utils import VariationType

def analyze_calls():
    num_cds, num_gene = 0, 0
    for receptor in gpcrdb.get_filtered_receptor_list("receptors.json"):
        entry_name = receptor['entry_name']
        receptor_class = receptor['receptor_class']
        dpath = os.path.join(receptor_class, entry_name)
        print(entry_name)

        calls_gene = set()
        with open(os.path.join(dpath, VCF_JPN_GENE_FILENAME)) as f:
            for l in f.readlines():
                cols = l.split('\t')
                chromosome = cols[0]
                position = int(cols[1])
                calls_gene.add((chromosome, position))

        calls_cds = set()
        with open(os.path.join(dpath, VCF_JPN_CDS_FILENAME)) as f:
            for l in f.readlines():
                cols = l.split('\t')
                chromosome = cols[0]
                position = int(cols[1])
                calls_cds.add((chromosome, position))

        assert(calls_cds.issubset(calls_gene))
        num_cds += len(calls_cds)
        num_gene += len(calls_gene)
    print("Total", num_cds, num_gene)

    fig, ax = plt.subplots(1, 1, figsize=(4, 2), dpi=300)

    ax.barh(0, num_cds, color='tab:orange')
    coding_text = "Coding region\n{:,} calls ({:.1f}%)".format(num_cds, num_cds / num_gene * 100)
    ax.text(num_cds / 2, 0.45, coding_text, ha='center', va='bottom', size=12)

    ax.barh(0, num_gene - num_cds, left=num_cds, color='tab:gray', alpha=0.6)
    non_coding_text = "Non-coding region\n{:,} calls ({:.1f}%)".format(num_gene - num_cds, (num_gene - num_cds) / num_gene * 100)
    ax.text((num_gene - num_cds) / 2, 0, non_coding_text, ha='center', va='center', size=12)

    ax.set_xlim(0, num_gene)
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig("./figures/S1a_calls.pdf")

def analyze_var_type():
    num_missense, num_silent, num_nonsense = 0, 0, 0
    for receptor in gpcrdb.get_filtered_receptor_list("receptors.json"):
        entry_name = receptor['entry_name']
        receptor_class = receptor['receptor_class']
        dpath = os.path.join(receptor_class, entry_name)
        with open(os.path.join(dpath, CSV_JPN_CDS_FILENAME)) as f:
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
    

    total = num_missense + num_silent + num_nonsense
    print(total, num_missense, num_silent, num_nonsense)

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
    fig.savefig("./figures/1a_var_type.pdf")

def analyze_var_seg():
    seg_missense, seg_silent, seg_nonsense = {}, {}, {}
    for receptor in gpcrdb.get_filtered_receptor_list("receptors.json"):
        entry_name = receptor['entry_name']
        receptor_class = receptor['receptor_class']
        dpath = os.path.join(receptor_class, entry_name)
        with open(os.path.join(dpath, CSV_JPN_CDS_FILENAME)) as f:
            for l in f.readlines():
                cols = l.split('\t')
                if l.startswith('#'):
                    header = cols
                    continue
                var_type = cols[header.index('#Var_Type')]
                seg = cols[header.index('Seg')]
                if seg == 'None':
                    print(dpath, l)
                if var_type == VariationType.MISSENSE.name:
                    seg_missense[seg] = seg_missense.get(seg, 0) + 1
                elif var_type == VariationType.SILENT.name:
                    seg_silent[seg] = seg_silent.get(seg, 0) + 1
                elif var_type == VariationType.NONSENSE.name:
                    seg_nonsense[seg] = seg_nonsense.get(seg, 0) + 1

    xlabels = ["N-term", "TM1", "ICL1", "TM2", "ECL1", "TM3", "ICL2", "TM4", "ECL2", "TM5", "ICL3", "TM6", "ECL3", "TM7", "H8", "ICL4", "C-term", "None"]
    cmap = plt.get_cmap('rainbow', len(xlabels) - 1)
    colors = [cmap(i) for i in range(len(xlabels) - 1)] + ['tab:gray']

    total_missense = sum([v for v in seg_missense.values()])
    total_silent = sum([v for v in seg_silent.values()])
    total_nonsense = sum([v for v in seg_nonsense.values()])

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
    ax.set_xlabel("SNVs [%]")
    fig.tight_layout()
    fig.savefig("./figures/1b_var_seg_percent.pdf")

if __name__ == '__main__':
    analyze_calls()
    analyze_var_type()
    analyze_var_seg()