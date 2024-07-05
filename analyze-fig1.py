#!/usr/bin/env python3
from config import *
import vcf
import gpcrdb
import matplotlib
matplotlib.use('Agg')
matplotlib.rc('pdf', fonttype=42)
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = "Arial"
import numpy as np
from utils import VariationType, Segment
import ensembl
import json

def analyze_calls(filename):
    num_cds, num_gene = 0, 0
    num_missense, num_silent, num_nonsense = 0, 0, 0
    for receptor in gpcrdb.get_filtered_receptor_list():
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

    fig, axes = plt.subplots(1, 2, figsize=(8, 2), dpi=300, sharey=True)

    ax = axes[0]
    num_gene_only = num_gene - num_cds
    ax.barh(0, num_gene_only, color='tab:gray', alpha=0.6)
    non_coding_text = "Non-coding region\n{:,} calls ({:.1f}%)".format(num_gene_only, num_gene_only / num_gene * 100)
    ax.text(num_gene_only / 2, 0, non_coding_text, ha='center', va='center')

    ax.barh(0, num_cds, left=num_gene_only, color='tab:orange')
    coding_text = "Coding region\n{:,} calls\n({:.1f}%)".format(num_cds, num_cds / num_gene * 100)
    ax.text(num_gene + num_cds / 2, -0.45, coding_text, ha='center', va='top')

    ax.set_xlim(0, num_gene)
    ax.set_axis_off()

    ax = axes[1]
    total = num_missense + num_silent + num_nonsense

    ax.barh(0, num_missense, color='tab:orange')
    missense_text = "Missense\n{:,} SNVs\n({:.1f}%)".format(num_missense, num_missense / total * 100)
    ax.text(num_missense / 2, 0, missense_text, ha='center', va='center')

    ax.barh(0, num_silent, left=num_missense, color='tab:gray', alpha=0.6)
    silent_text = "Silent\n{:,} SNVs\n({:.1f}%)".format(num_silent, num_silent / total * 100)
    ax.text(num_missense + num_silent / 2, 0, silent_text, ha='center', va='center')

    ax.barh(0, num_nonsense, left=num_missense + num_silent, color='tab:gray')
    nonsense_text = "Nonsense\n{:,} SNVs\n({:.1f}%)".format(num_nonsense, num_nonsense / total * 100)
    ax.text(num_missense + num_silent + num_nonsense / 2, -0.45, nonsense_text, ha='center', va='top')

    ax.set_xlim(0, total)
    ax.set_axis_off()

    fig.tight_layout()
    fig.savefig(filename)

def analyze_var_seg(filename):
    seg_missense, seg_silent, seg_nonsense = {}, {}, {}
    for receptor in gpcrdb.get_filtered_receptor_list():
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
    fig.savefig(filename)

def analyze_segment_ratio(filename_A, filename_B):
    residue_counts = {}
    for receptor in gpcrdb.get_filtered_receptor_list():
        counts = {seg: 0 for seg in Segment}
        with open(receptor.alignment_path) as f:
            for l in f:
                if l.startswith('#'):
                    continue
                cols = l.split(',')
                seg = Segment.value_of(cols[1])
                counts[seg] = counts.get(seg, 0) + 1

        d = residue_counts.get(receptor.receptor_class, {seg: [] for seg in Segment})
        for seg in Segment:
            d[seg].append(counts[seg])
        residue_counts[receptor.receptor_class] = d

    receptor_classes = sorted(residue_counts.keys(), reverse=True)
    
    fig, ax = plt.subplots(1, 1, figsize=(3, 4), dpi=300)

    ticklabels = []
    for y, family in enumerate(receptor_classes):
        left = 0
        total = sum([sum(residue_counts[family][seg]) for seg in Segment])
        num = max([len(residue_counts[family][seg]) for seg in Segment])
        family_label = "Family " + family.split(' ')[1] if family.startswith('Class') else "Others"
        ticklabels.append(family_label + "\nn = {}".format(num))

        for seg in Segment:
            width = sum(residue_counts[family][seg]) / total * 100
            ax.barh(y, width, height=0.7, left=left, color=seg.color, edgecolor='k', lw=0.2)
            left += width

    left = 0
    total = sum([sum([sum(residue_counts[family][seg]) for seg in Segment]) for family in receptor_classes])
    ticklabels.append("All")

    for seg in Segment:
        y = len(receptor_classes)
        width = sum([sum(residue_counts[family][seg]) for family in receptor_classes]) / total * 100
        ax.barh(y, width, height=0.7, left=left, color=seg.color, edgecolor='k', lw=0.2)

        if seg in Segment.TMs() or seg in Segment.terms():
            if seg == Segment.TM1 or seg in Segment.terms():
                ax.text(left + width / 2, y + 0.4, seg.value, ha='center', va='bottom', size=5)
            else:
                ax.text(left + width / 2, y + 0.4, seg.value[2], ha='center', va='bottom', size=5)
        left += width

    ax.set_xlabel("Segment fractions [%]")
    ax.set_yticks(range(len(receptor_classes) + 1))
    ax.set_yticklabels(ticklabels, size=8)
    
    fig.tight_layout()
    fig.savefig(filename_A)
    plt.close(fig)

    fig, axes = plt.subplots(4, 4, figsize=(8, 7), dpi=300, sharey=True)
    segments = [seg for seg in Segment if seg != Segment.ICL4]

    for y in range(4):
        for x in range(4):
            ax = axes[y][x]
            seg = segments[4 * y + x]

            values = sum([by_family[seg] for by_family in residue_counts.values()], [])
            bins = 20
            edges = (20, 60)
            if min(values) < edges[0] or edges[1] < max(values):
                edges = None

            data_A = residue_counts['Class A (Rhodopsin)'][seg]
            data_nonA = sum([residue_counts[family][seg] for family in residue_counts.keys() if family != 'Class A (Rhodopsin)'], [])
            ax.hist([data_A, data_nonA], bins=bins, range=edges, color=[seg.color, 'lightgray'], stacked=True)
            median = np.median(values)
            ax.axvline(median, color='k', lw=1)
            ax2 = ax.twiny()
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks([median])
            ax2.set_xticklabels([f"{median:.1f}"])
            ax2.tick_params(length=0, pad=0)

            ax.set_title(seg.value)
            ax.set_ylim(10**-0.5, 10**2.8)
            ax.set_yscale('log')
            if x == 0:
                ax.set_ylabel("Number of GPCRs")
            if y == 3:
                ax.set_xlabel("Number of residues [aa]")
    fig.tight_layout()
    fig.savefig(filename_B)

def analyze_gene_stats(filename):
    gene_lengths_A, gene_lengths_nonA = [], []
    transcipt_lengths_A, transcipt_lengths_nonA = [], []
    exon_numbers_A, exon_numbers_nonA = [], []
    translation_lengths_A, translation_lengths_nonA = [], []
    labels = ["Family A", "Others"]
    colors = ["tab:gray", "lightgray"]

    for receptor in gpcrdb.get_filtered_receptor_list():
        with open(receptor.ensembl_path) as f:
            j = json.load(f)

            start = j['start']
            end = j['end']
            if receptor.receptor_class == 'Class A (Rhodopsin)':
                gene_lengths_A.append(end - start + 1)
            else:
                gene_lengths_nonA.append(end - start + 1)

            for transcipt in j['Transcript']:
                if transcipt['is_canonical'] == 1:
                    exons = transcipt['Exon']

                    if receptor.receptor_class == 'Class A (Rhodopsin)':
                        transcipt_lengths_A.append(sum([exon['end'] - exon['start'] + 1 for exon in exons]))
                        exon_numbers_A.append(len(exons))
                        translation_lengths_A.append(transcipt['Translation']['length'])
                    else:
                        transcipt_lengths_nonA.append(sum([exon['end'] - exon['start'] + 1 for exon in exons]))
                        exon_numbers_nonA.append(len(exons))
                        translation_lengths_nonA.append(transcipt['Translation']['length'])

    fig, axes = plt.subplots(2, 2, figsize=(5, 4), dpi=300)
    
    ax = axes[0][0]
    ax.hist([np.log10(gene_lengths_A), np.log10(gene_lengths_nonA)], bins=40, color=colors, label=labels, stacked=True)
    ax.set_xlabel("Gene length [nt]")
    xrange = np.arange(3, 7)
    ax.set_xticks(xrange)
    ax.set_xticklabels([f"10$^{{{v}}}$" for v in xrange])
    median = np.median(np.concatenate([gene_lengths_A, gene_lengths_nonA]))
    ax.axvline(np.log10(median), color='tab:orange')
    ax2 = ax.twiny()
    ax2.set_xticks([np.log10(median)])
    ax2.set_xticklabels([f"{median:.1f}"])
    ax2.set_xlim(ax.get_xlim())
    ax2.tick_params(colors='tab:orange', length=0, pad=0)
    ax.set_ylabel("Number of GPCRs")

    ax = axes[1][0]
    ax.hist([np.log10(transcipt_lengths_A), np.log10(transcipt_lengths_nonA)], bins=40, color=colors, label=labels, stacked=True)
    ax.set_xlabel("Transcript length [nt]")
    xrange = np.arange(3, 4.5, 0.5)
    ax.set_xticks(xrange)
    ax.set_xticklabels([f"10$^{{{v}}}$" for v in xrange])
    median = np.median(np.concatenate([transcipt_lengths_A, transcipt_lengths_nonA]))
    ax.axvline(np.log10(median), color='tab:orange')
    ax2 = ax.twiny()
    ax2.set_xticks([np.log10(median)])
    ax2.set_xticklabels([f"{median:.1f}"])
    ax2.set_xlim(ax.get_xlim())
    ax2.tick_params(colors='tab:orange', length=0, pad=0)
    ax.set_ylabel("Number of GPCRs")

    ax = axes[0][1]
    ax.hist([exon_numbers_A, exon_numbers_nonA], bins=40, color=colors, label=labels, stacked=True)
    ax.set_xlabel("Number of exons")
    median = np.median(np.concatenate([exon_numbers_A, exon_numbers_nonA]))
    ax.axvline(median, color='tab:orange')
    ax2 = ax.twiny()
    ax2.set_xticks([median])
    ax2.set_xticklabels([f"{median:.1f}"])
    ax2.set_xlim(ax.get_xlim())
    ax2.tick_params(colors='tab:orange', length=0, pad=0)
    ax.set_yscale('log')
    ax.legend()
    
    ax = axes[1][1]
    ax.hist([np.log10(translation_lengths_A), np.log10(translation_lengths_nonA)], bins=40, color=colors, label=labels, stacked=True)
    ax.set_xlabel("Translation length [aa]")
    xrange = np.arange(2.5, 4, 0.5)
    ax.set_xticks(xrange)
    ax.set_xticklabels([f"10$^{{{v}}}$" for v in xrange])
    median = np.median(np.concatenate([translation_lengths_A, translation_lengths_nonA]))
    ax.axvline(np.log10(median), color='tab:orange')
    ax2 = ax.twiny()
    ax2.set_xticks([np.log10(median)])
    ax2.set_xticklabels([f"{median:.1f}"])
    ax2.set_xlim(ax.get_xlim())
    ax2.tick_params(colors='tab:orange', length=0, pad=0)
    ax.set_yscale('log')

    fig.tight_layout()
    fig.savefig(filename)

if __name__ == '__main__':
    analyze_calls("./figures/1ab_variations.pdf")
    analyze_var_seg("./figures/1d_var_seg_percent.pdf")
    analyze_calls("./figures/1ab_variations.pdf")
    # analyze_gene_stats("./figures/S1a_stats.pdf")
    # analyze_segment_ratio("./figures/S1b_segment_ratio.pdf", "./figures/S1c_residue_count.pdf")