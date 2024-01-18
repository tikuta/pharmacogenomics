#!/usr/bin/env python3
from config import *
import vcf
import gpcrdb
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = "Arial"
import numpy as np
from utils import VariationType, Segment
import ensembl
import json

def analyze_calls():
    num_cds, num_gene = 0, 0
    for receptor in gpcrdb.get_filtered_receptor_list("receptors.json"):
        calls_gene = set()
        with open(receptor.japan_gene_vcf_path) as f:
            for l in f.readlines():
                try:
                    var = vcf.VariationEntry.load_from_54KJPN(l)
                    calls_gene.add((var.chromosome, var.position))
                except (vcf.NotPassedError, vcf.BlankLineError):
                    continue

        calls_cds = set()
        with open(receptor.japan_cds_vcf_path) as f:
            for l in f.readlines():
                try:
                    var = vcf.VariationEntry.load_from_54KJPN(l)
                    calls_cds.add((var.chromosome, var.position))
                except (vcf.NotPassedError, vcf.BlankLineError):
                    continue

        assert(calls_cds.issubset(calls_gene))
        num_cds += len(calls_cds)
        num_gene += len(calls_gene)

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
    fig.savefig("./figures/S1a_calls.pdf")

def analyze_var_type():
    num_missense, num_silent, num_nonsense = 0, 0, 0
    for receptor in gpcrdb.get_filtered_receptor_list("receptors.json"):
        with open(receptor.japan_cds_csv_path) as f:
            for l in f.readlines():
                try:
                    anno = ensembl.Annotation.from_csv_line(l)
                    if anno.var_type == VariationType.MISSENSE:
                        num_missense += 1
                    elif anno.var_type == VariationType.SILENT:
                        num_silent += 1
                    elif anno.var_type == VariationType.NONSENSE:
                        num_nonsense += 1
                except ensembl.BlankLineError:
                    continue

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
    fig.savefig("./figures/1a_var_type.pdf")

def analyze_var_seg():
    seg_missense, seg_silent, seg_nonsense = {}, {}, {}
    for receptor in gpcrdb.get_filtered_receptor_list("receptors.json"):
        with open(receptor.japan_cds_csv_path) as f:
            for l in f.readlines():
                try:
                    anno = ensembl.Annotation.from_csv_line(l)
                    if anno.segment == None:
                        print(receptor.entry_name, anno.snv)
                        continue

                    if anno.var_type == VariationType.MISSENSE:
                        seg_missense[anno.segment.value] = seg_missense.get(anno.segment.value, 0) + 1
                    elif anno.var_type == VariationType.SILENT:
                        seg_silent[anno.segment.value] = seg_silent.get(anno.segment.value, 0) + 1
                    elif anno.var_type == VariationType.NONSENSE:
                        seg_nonsense[anno.segment.value] = seg_nonsense.get(anno.segment.value, 0) + 1
                except ensembl.BlankLineError:
                    continue

    total_missense = sum([v for v in seg_missense.values()])
    total_silent = sum([v for v in seg_silent.values()])
    total_nonsense = sum([v for v in seg_nonsense.values()])

    fig, ax = plt.subplots(1, 1, figsize=(4, 2), dpi=300)

    left = 0
    for seg in Segment:
        width = seg_missense.get(seg.value, 0) / total_missense * 100
        ax.barh(2, width, height=0.7, left=left, color=seg.color, edgecolor='k', lw=0.2)
        
        label = seg.value
        if label.startswith('TM') or label.endswith('-term'):
            ax.text(left + width / 2, 2.4, label, ha='center', va='bottom', size=6)
        left += width

    left = 0
    for seg in Segment:
        width = seg_silent.get(seg.value, 0) / total_silent * 100
        ax.barh(1, width, height=0.7, left=left, color=seg.color, edgecolor='k', lw=0.2)
        left += width

    left = 0
    for seg in Segment:
        width = seg_nonsense.get(seg.value, 0) / total_nonsense * 100
        ax.barh(0, width, height=0.7, left=left, color=seg.color, edgecolor='k', lw=0.2)
        left += width
    num_loss_of_function = sum([seg_nonsense.get(seg.value, 0) for seg in Segment if seg not in (Segment.H8, Segment.Cterm)])
    ratio_loss_of_function = num_loss_of_function / total_nonsense * 100
    ax.plot([0, ratio_loss_of_function], [-0.45, -0.45], lw=1, color='tab:gray')
    ax.text(ratio_loss_of_function / 2, -0.5, "{:.1f}%".format(ratio_loss_of_function), ha='center', va='top', size=6, color='tab:gray')

    ax.set_ylim(-0.75, 2.7)
    ax.set_yticks(range(3))
    ax.set_yticklabels(["Nonsense", "Silent", "Missense"])
    ax.set_xlabel("SNVs [%]")
    fig.tight_layout()
    fig.savefig("./figures/1b_var_seg_percent.pdf")

def analyze_segment_ratio():
    residue_counts = {seg: [] for seg in Segment}
    for receptor in gpcrdb.get_filtered_receptor_list("receptors.json"):
        counts = {seg: 0 for seg in Segment}
        with open(receptor.alignment_path) as f:
            for l in f:
                if l.startswith('#'):
                    continue
                cols = l.split(',')
                seg = Segment.value_of(cols[1])
                counts[seg] = counts.get(seg, 0) + 1
        for seg in residue_counts.keys():
            residue_counts[seg].append(counts[seg])
    
    fig, axes = plt.subplots(8, 2, figsize=(10, 8), dpi=300)
    helices = [seg for seg in Segment if seg.value.startswith('TM') or seg == Segment.H8]
    terms = [Segment.Nterm, Segment.Cterm]
    loops = [Segment.ICL1, Segment.ICL2, Segment.ICL3, Segment.ECL1, Segment.ECL2, Segment.ECL3]

    medianprops = {'color': 'black'}
    meanprops = {'markerfacecolor': 'black', 'markeredgewidth': 0, 'marker': '.', 'markersize': 10}
    flierprops = {'markeredgecolor': 'black', 'markeredgewidth': 1, 'marker': 'x', 'markersize': 3}

    for ax, seg in zip(axes.T[0], helices):
        boxprops = {'color': 'black', 'facecolor': seg.color}

        ax.boxplot(residue_counts[seg], patch_artist=True, labels=[seg.value], vert=False, showmeans=True,
                    medianprops=medianprops, flierprops=flierprops, boxprops=boxprops, meanprops=meanprops)
        
        min, median, mean, max = np.min(residue_counts[seg]), np.median(residue_counts[seg]), np.mean(residue_counts[seg]), np.max(residue_counts[seg])
        ax.text(min, 1, "{:.0f} ".format(min), ha='right', va='center')
        ax.text(max, 1, " {:.0f}".format(max), ha='left', va='center')
        ax.text(median, 1.1, "{:.1f}".format(median), ha='center', va='bottom')
        ax.text(mean, 0.85, "{:.1f}".format(mean), ha='center', va='top')
        ax.set_ylim(0.7, 1.3)
        if seg == Segment.H8:
            ax.set_xlim(-2, 30)
        else:
            ax.set_xlim(18, 62)
    
    right_limits = {
        Segment.Nterm: 7000, Segment.Cterm: 2000, 
        Segment.ICL1: 17, Segment.ICL2: 27, Segment.ICL3: 240,
        Segment.ECL1: 67, Segment.ECL2: 200, Segment.ECL3: 33
    }
    for ax, seg in zip(axes.T[1], terms + loops):
        boxprops = {'color': 'black', 'facecolor': seg.color}

        ax.boxplot(residue_counts[seg], patch_artist=True, labels=[seg.value], vert=False, showmeans=True,
                    medianprops=medianprops, flierprops=flierprops, boxprops=boxprops, meanprops=meanprops)
        
        min, median, mean, max = np.min(residue_counts[seg]), np.median(residue_counts[seg]), np.mean(residue_counts[seg]), np.max(residue_counts[seg])
        ax.text(min, 1, "{:.0f} ".format(min), ha='right', va='center')
        ax.text(max, 1, " {:.0f}".format(max), ha='left', va='center')
        ax.text(median, 1.1, "{:.1f}".format(median), ha='center', va='bottom')
        ax.text(mean, 0.85, "{:.1f}".format(mean), ha='center', va='top')
        ax.set_ylim(0.7, 1.3)
        ax.set_xlim(right=right_limits[seg])
    axes[-1][0].set_xlabel("Number of residues")
    axes[-1][1].set_xlabel("Number of residues")
    fig.tight_layout()
    fig.savefig("./figures/S1d_residue_count.pdf")

def analyze_gene_stats():
    gene_lengths = []
    transcipt_lengths = []
    exon_numbers = []
    translation_lengths = []

    for receptor in gpcrdb.get_filtered_receptor_list("receptors.json"):
        with open(receptor.ensembl_path) as f:
            j = json.load(f)

            start = j['start']
            end = j['end']
            gene_lengths.append(end - start + 1)

            for transcipt in j['Transcript']:
                if transcipt['is_canonical'] == 1:
                    exons = transcipt['Exon']
                    transcipt_lengths.append(sum([exon['end'] - exon['start'] + 1 for exon in exons]))
                    exon_numbers.append(len(exons))
                    translation_lengths.append(transcipt['Translation']['length'])

    boxprops = {'facecolor': 'tab:gray', 'color': 'tab:gray'}
    whisker_cap_props = {'color': 'tab:gray'}
    medianprops = {'color': 'tab:orange'}
    meanprops = {'markerfacecolor': 'tab:orange', 'markeredgecolor': 'tab:orange', 'marker': '.'}
    flierprops = {'markeredgecolor': 'tab:gray', 'markeredgewidth': 1, 'marker': 'x', 'markersize': 3}

    fig, axes = plt.subplots(4, 1, figsize=(7, 5), dpi=300)

    axes[0].boxplot(gene_lengths, vert=False, patch_artist=True, showmeans=True, 
               medianprops=medianprops, meanprops=meanprops, boxprops=boxprops, 
               whiskerprops=whisker_cap_props, capprops=whisker_cap_props, flierprops=flierprops)
    min, median, mean, max = np.min(gene_lengths), np.median(gene_lengths), np.mean(gene_lengths), np.max(gene_lengths)
    axes[0].text(min, 1, "{:.0f} ".format(min), ha='right', va='center', color='tab:gray')
    axes[0].text(max, 1, " {:.0f}".format(max), ha='left', va='center', color='tab:gray')
    axes[0].text(median, 1.1, "{:.1f}".format(median), ha='center', va='bottom', color='tab:orange')
    axes[0].text(mean, 0.85, "{:.1f}".format(mean), ha='center', va='top', color='tab:orange')
    axes[0].set_xlabel("Gene length [nucleotides]")
    axes[0].set_xscale('log')
    axes[0].set_xlim(10**2.5, 10**6.5)
    axes[0].set_yticks([])

    axes[1].boxplot(transcipt_lengths, vert=False, patch_artist=True, showmeans=True, 
               medianprops=medianprops, meanprops=meanprops, boxprops=boxprops, 
               whiskerprops=whisker_cap_props, capprops=whisker_cap_props, flierprops=flierprops)
    min, median, mean, max = np.min(transcipt_lengths), np.median(transcipt_lengths), np.mean(transcipt_lengths), np.max(transcipt_lengths)
    axes[1].text(min, 1, "{:.0f} ".format(min), ha='right', va='center', color='tab:gray')
    axes[1].text(max, 1, " {:.0f}".format(max), ha='left', va='center', color='tab:gray')
    axes[1].text(median, 1.1, "{:.1f}".format(median), ha='center', va='bottom', color='tab:orange')
    axes[1].text(mean, 0.85, "{:.1f}".format(mean), ha='center', va='top', color='tab:orange')
    axes[1].set_xlabel("Canonical transcript length [nucleotides]")
    axes[1].set_xlim(-1000, 23000)
    axes[1].set_yticks([])

    axes[2].boxplot(exon_numbers, vert=False, patch_artist=True, showmeans=True, 
               medianprops=medianprops, meanprops=meanprops, boxprops=boxprops, 
               whiskerprops=whisker_cap_props, capprops=whisker_cap_props, flierprops=flierprops)
    min, median, mean, max = np.min(exon_numbers), np.median(exon_numbers), np.mean(exon_numbers), np.max(exon_numbers)
    axes[2].text(min, 1, "{:.0f} ".format(min), ha='right', va='center', color='tab:gray')
    axes[2].text(max, 1, " {:.0f}".format(max), ha='left', va='center', color='tab:gray')
    axes[2].text(median, 1.1, "{:.1f}".format(median), ha='center', va='bottom', color='tab:orange')
    axes[2].text(mean, 0.85, "{:.1f}".format(mean), ha='center', va='top', color='tab:orange')
    axes[2].set_xlabel("Number of exons")
    axes[2].set_xlim(-5, 95)
    axes[2].set_yticks([])

    axes[3].boxplot(translation_lengths, vert=False, patch_artist=True, showmeans=True, 
               medianprops=medianprops, meanprops=meanprops, boxprops=boxprops, 
               whiskerprops=whisker_cap_props, capprops=whisker_cap_props, flierprops=flierprops)
    min, median, mean, max = np.min(translation_lengths), np.median(translation_lengths), np.mean(translation_lengths), np.max(translation_lengths)
    axes[3].text(min, 1, "{:.0f} ".format(min), ha='right', va='center', color='tab:gray')
    axes[3].text(max, 1, " {:.0f}".format(max), ha='left', va='center', color='tab:gray')
    axes[3].text(median, 1.1, "{:.1f}".format(median), ha='center', va='bottom', color='tab:orange')
    axes[3].text(mean, 0.85, "{:.1f}".format(mean), ha='center', va='top', color='tab:orange')   
    axes[3].set_xlabel("Canonical translation length [amino acids]")
    axes[3].set_xlim(-200, 7000)
    axes[3].set_yticks([])

    fig.tight_layout()
    fig.savefig("./figures/S1b_stats.pdf")

if __name__ == '__main__':
    analyze_calls()
    analyze_var_type()
    analyze_var_seg()
    analyze_gene_stats()
    analyze_segment_ratio()