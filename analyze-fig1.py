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

    ax.set_ylim(-0.5, 2.7)
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
    
    fig, ax = plt.subplots(1, 1, figsize=(8, 5), dpi=300)
    segments = [seg for seg in Segment if seg not in (Segment.NONE, Segment.FailedToGuess)]
    labels = ["{}\n({})".format(seg.value, np.median(residue_counts[seg])) for seg in segments]
    bplot = ax.boxplot([residue_counts[seg] for seg in segments], sym='', patch_artist=True, labels=labels)
    for box, median, seg in zip(bplot['boxes'], bplot['medians'], segments):
        box.set_facecolor(seg.color)
        median.set_color('black')
    ax.set_ylim(bottom=-5)
    ax.set_ylabel("Number of residues")
    ax.set_xlabel("Segment\n(median)")
    fig.tight_layout()
    fig.savefig("./figures/S1c_residue_count.pdf")

if __name__ == '__main__':
    analyze_calls()
    analyze_var_type()
    analyze_var_seg()
    analyze_segment_ratio()