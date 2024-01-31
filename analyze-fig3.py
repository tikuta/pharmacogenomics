#!/usr/bin/env python3
import gpcrdb
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = "Arial"
import ensembl
from utils import VariationType, Segment

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

if __name__ == '__main__':
    analyze_positions()