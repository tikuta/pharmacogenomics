#!/usr/bin/env python3
import gpcrdb
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = "Arial"
import ensembl
from utils import VariationType, Segment
import json
import config

def analyze_positions():
    assigned = {}
    found = {}
    for receptor in gpcrdb.get_filtered_receptor_list("receptors.json"):
        if receptor.receptor_class != 'Class A (Rhodopsin)':
            continue

        # Generic number is unique in a receptor
        with open(receptor.alignment_path) as f:
            for l in f:
                if l.startswith('#'):
                    continue
                cols = l.split(',')
                generic_number = cols[2]
                if generic_number == 'None':
                    continue
                latter = generic_number.split('x')[-1]
                if len(latter) > 2: # Insertion
                    continue
                structure_based_number = generic_number.split('.')[0] + 'x' + latter
                assigned[structure_based_number] = assigned.get(structure_based_number, 0) + 1
        
        # Same generic number could appear in a CSV file 
        missense = set()
        with open(receptor.japan_cds_csv_path) as f:
            for l in f.readlines():
                try:
                    anno = ensembl.Annotation.from_csv_line(l)
                except ensembl.BlankLineError:
                    continue

                if anno.var_type != VariationType.MISSENSE:
                    continue

                if not anno.generic_number:
                    continue
                
                latter = anno.generic_number.split('x')[-1]
                if len(latter) > 2: # Insertion
                    continue
                structure_based_number = anno.generic_number.split('.')[0] + 'x' + latter
                missense.add(structure_based_number)
        
        for gn in missense:
            found[gn] = found.get(gn, 0) + 1

    f = lambda gn: (int(gn.split('x')[0]) if int(gn.split('x')[0]) > 10 else int(gn.split('x')[0]) * 10, int(gn.split('x')[1]))
    generic_numbers = sorted(assigned.keys(), key=f)

    fig, ax = plt.subplots(1, 1, figsize=(4, 10), dpi=300)
    ax.invert_yaxis()
    
    for y, gn in enumerate(generic_numbers):
        color = Segment.generic_number_of(gn).color
        ax.barh(y, assigned[gn], color='tab:gray', height=1, zorder=-200)
        ax.barh(y, found.get(gn, 0), color=color, height=1, zorder=-100)

    ax.set_ylim(len(assigned), -1)
    ax.set_yticks([y for y, gn in enumerate(generic_numbers) if gn.endswith('x50')])
    ax.set_yticklabels([gn for gn in generic_numbers if gn.endswith('x50')])
    ax.set_ylabel("Structure-based generic number")
    num_receptors = max(assigned.values())
    ax.set_xlim(0, num_receptors)
    ax.set_xlabel("Number of family A GPCRs")
    ax.set_xticks([0, 50, 100, 150, 200, 250, num_receptors])
    ax.set_xticklabels([0, 50, 100, 150, 200, 250, num_receptors])
    fig.tight_layout()
    fig.savefig("./figures/3a_positions.pdf")

def analyze_arginine_3x50():
    roi = "3x50"

    aa_stats = {}
    codon_stats = {}
    for receptor in gpcrdb.get_filtered_receptor_list("receptors.json"):
        if receptor.receptor_class != 'Class A (Rhodopsin)':
            continue

        with open(receptor.alignment_path) as f:
            for l in f:
                if l.startswith('#'):
                    continue

                cols = l.strip().split(',')
                residue, generic_number = cols[0], cols[2]
                structure_based_number = generic_number.split('.')[0] + 'x' + generic_number.split('x')[-1]

                if structure_based_number == roi:
                    aa = residue[0]
                    aa_stats[aa] = aa_stats.get(aa, 0) + 1

                    res_num = int(residue[1:])

                    if aa == 'R':
                        with open(receptor.cds_path) as j:
                            codon = json.load(j)['canonical_cds'][res_num * 3 - 3: res_num * 3]
                            codon_stats[codon] = codon_stats.get(codon, 0) + 1
    
    fig, ax = plt.subplots(1, 1, figsize=(6, 2), dpi=300)

    left = 0
    for aa in sorted(aa_stats.keys(), key=lambda aa:aa_stats[aa], reverse=True):
        delta = aa_stats[aa]
        ax.barh(1, delta, left=left, color=config.AA2COLOR[aa], linewidth=0.5, edgecolor='k')
        if aa == "R":
            ax.text(left + delta / 2, 1, "Arg\n({})".format(delta), ha='center', va='center', size=6)
        left += delta
    total = left
    
    left = 0
    for codon in sorted(codon_stats.keys(), key=lambda codon: ('CG' not in codon, -codon_stats[codon])):
        delta = codon_stats[codon]
        ax.barh(0, delta, left=left, linewidth=0.5, edgecolor='k')
        ax.text(left + delta / 2, 0, "{}\n({})".format(codon, delta), ha='center', va='center', size=6)
        left += delta
    cpg = sum([v for k, v in codon_stats.items() if 'CG' in k])

    ax.set_xlim(0, total)
    ax.set_xticks([0, 50, 100, 150, 200, 250, total])
    ax.set_xticklabels([0, 50, 100, 150, 200, 250, total])
    ax.set_xlabel("Number of family A GPCRs")

    ax.set_yticks([0, 1])
    ax.set_yticklabels(['Codon', 'Amino acid'])
    ax.set_ylabel("3x50")
    
    fig.tight_layout()
    fig.savefig("./figures/S3a_arginine_3x50.pdf")

if __name__ == '__main__':
    analyze_positions()
    analyze_arginine_3x50()