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
    for receptor in gpcrdb.get_filtered_receptor_list():
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
    fig.savefig("./figures/1b_var_seg_percent.pdf")

def analyze_segment_ratio():
    residue_counts = {seg: [] for seg in Segment}
    for receptor in gpcrdb.get_filtered_receptor_list():
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
    
    fig, axes = plt.subplots(4, 4, figsize=(8, 7), dpi=300, sharey=True)
    segments = [seg for seg in Segment if seg != Segment.ICL4]

    for y in range(4):
        for x in range(4):
            ax = axes[y][x]
            seg = segments[4 * y + x]

            values = residue_counts[seg]
            bins = 20
            edges = (20, 60)
            if np.min(values) < edges[0] or edges[1] < np.max(values):
                edges = None

            ax.hist(values, bins=bins, range=edges, color=seg.color)
            median = np.median(values)
            ax.axvline(median, color='k', lw=1)

            title = seg.value
            if seg == Segment.Nterm:
                title += f" (median = {median})"
            else:
                title += f" ({median})"
            ax.set_title(title)
            ax.set_yscale('log')
            if x == 0:
                ax.set_ylabel("Number of GPCRs")
            if y == 3:
                ax.set_xlabel("Number of residues\n[amino acids]")
    fig.tight_layout()
    fig.savefig("./figures/S1d_residue_count.pdf")

def analyze_gene_stats():
    gene_lengths = []
    transcipt_lengths = []
    exon_numbers = []
    translation_lengths = []

    for receptor in gpcrdb.get_filtered_receptor_list():
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

    fig, axes = plt.subplots(2, 2, figsize=(7, 5), dpi=300)
    
    nbins, _, _ = axes[0][0].hist(np.log10(gene_lengths), bins=40, color='tab:gray')
    axes[0][0].set_xlabel("Gene length [nucleotides]")
    xrange = np.arange(3, 7)
    axes[0][0].set_xticks(xrange)
    axes[0][0].set_xticklabels([f"10$^{{{v}}}$" for v in xrange])
    median = np.median(gene_lengths)
    axes[0][0].axvline(np.log10(median), color='tab:orange')
    axes[0][0].text(np.log10(median), np.max(nbins), f"{median:.0f} ", color='tab:orange', ha='right', va='center')
    axes[0][0].set_ylabel("Number of GPCRs")

    nbins, _, _ = axes[1][0].hist(np.log10(transcipt_lengths), bins=40, color='tab:gray')
    axes[1][0].set_xlabel("Canonical transcript length [nucleotides]")
    xrange = np.arange(3, 4.5, 0.25)
    axes[1][0].set_xticks(xrange)
    axes[1][0].set_xticklabels([f"10$^{{{v}}}$" for v in xrange])
    median = np.median(transcipt_lengths)
    axes[1][0].axvline(np.log10(median), color='tab:orange')
    axes[1][0].text(np.log10(median), np.max(nbins), f" {median:.0f}", color='tab:orange', ha='left', va='center')
    axes[1][0].set_ylabel("Number of GPCRs")

    nbins, _, _ = axes[0][1].hist(exon_numbers, bins=20, color='tab:gray')
    axes[0][1].set_xlabel("Number of exons")
    median = np.median(exon_numbers)
    axes[0][1].axvline(median, color='tab:orange')
    axes[0][1].text(np.log10(median), np.max(nbins), f"    {median:.0f}", color='tab:orange', ha='left', va='center')
    axes[0][1].set_yscale('log')
    
    nbins, _, _ = axes[1][1].hist(np.log10(translation_lengths), bins=20, color='tab:gray')
    axes[1][1].set_xlabel("Canonical translation length [amino acids]")
    xrange = np.arange(2.5, 4, 0.25)
    axes[1][1].set_xticks(xrange)
    axes[1][1].set_xticklabels([f"10$^{{{v}}}$" for v in xrange])
    median = np.median(translation_lengths)
    axes[1][1].axvline(np.log10(median), color='tab:orange')
    axes[1][1].text(np.log10(median), np.max(nbins), f" {median:.0f}", color='tab:orange', ha='left', va='center')
    axes[1][1].set_yscale('log')

    fig.tight_layout()
    fig.savefig("./figures/S1b_stats.pdf")

if __name__ == '__main__':
    # analyze_calls()
    # analyze_var_type()
    # analyze_var_seg()
    # analyze_gene_stats()
    analyze_segment_ratio()