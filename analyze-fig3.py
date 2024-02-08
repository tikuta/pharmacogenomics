#!/usr/bin/env python3
import gpcrdb
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = "Arial"
import ensembl
from utils import VariationType, Segment, GproteinCoupling
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

    generic_numbers = sorted(assigned.keys(), key=lambda gn:Segment.generic_number_of(gn).index)

    fig, ax = plt.subplots(1, 1, figsize=(4, 10), dpi=300)
    ax.invert_yaxis()
    
    for y, gn in enumerate(generic_numbers):
        color = Segment.generic_number_of(gn).color
        ax.barh(y, assigned[gn], color='lightgray', height=1, zorder=-200)
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

    ax.set_xlim(0, total)
    ax.set_xticks([0, 50, 100, 150, 200, 250, total])
    ax.set_xticklabels([0, 50, 100, 150, 200, 250, total])
    ax.set_xlabel("Number of family A GPCRs")

    ax.set_yticks([0, 1])
    ax.set_yticklabels(['Codon', 'Amino acid'])
    ax.set_ylabel("3x50")
    
    fig.tight_layout()
    fig.savefig("./figures/S3a_arginine_3x50.pdf")

def analyze_G_protein_contact_positions():
    # See the following reference for detail
    # https://doi.org/10.1038/s41467-022-34055-5
    common_residues = {"3x50", "3x53", "3x54", "34x50", "34x51", "34x55", "5x65", "5x68", "6x32", "6x33", "6x36", "6x37", "7x56", "8x47"}
    gs_residues = common_residues | {"3x54", "3x55", "34x51", "34x54", "34x55", "5x64", "5x68", "5x69", "5x71", "5x72", "5x74", "5x75", "5x77", "8x48"}
    gi_residues = common_residues | {"12x49", "2x40", "3x50", "3x53", "34x52", "34x55", "5x71", "6x32", "7x56", "8x47", "8x49"}
    gq_residues = common_residues | {"2x37", "2x39", "2x40", "3x49", "34x51", "34x53", "34x55", "34x56", "34x57", "4x38", "4x39", "6x30", "6x33", "8x48", "8x49"}

    roi = gs_residues | gi_residues | gq_residues
    gs_gi_gq = gs_residues & gi_residues & gq_residues
    gs = gs_residues - gi_residues - gq_residues
    gi = gi_residues - gs_residues - gq_residues
    gq = gq_residues - gs_residues - gi_residues
    gs_gi = gs_residues & gi_residues - gq_residues
    gi_gq = gi_residues & gq_residues - gs_residues
    gq_gs = gq_residues & gs_residues - gi_residues

    # CMYK coloring (C = Gs, M = Gi, Y = Gq, K = common)
    colormap = {}
    print("All", " ".join(sorted(list(roi), key=lambda gn: Segment.generic_number_of(gn).index)))
    print("Common (gray)", " ".join(sorted(list(gs_gi_gq), key=lambda gn: Segment.generic_number_of(gn).index)))
    colormap.update({gn: "gray50" for gn in gs_gi_gq})
    print("Gs (cyan)", " ".join(sorted(list(gs), key=lambda gn: Segment.generic_number_of(gn).index)))
    colormap.update({gn: "cyan" for gn in gs})
    print("Gi (magenta)", " ".join(sorted(list(gi), key=lambda gn: Segment.generic_number_of(gn).index)))
    colormap.update({gn: "magenta" for gn in gi})
    print("Gq (yellow)", " ".join(sorted(list(gq), key=lambda gn: Segment.generic_number_of(gn).index)))
    colormap.update({gn: "yellow" for gn in gq})
    print("Gs^Gi (blue)", " ".join(sorted(list(gs_gi), key=lambda gn: Segment.generic_number_of(gn).index)))
    colormap.update({gn: "blue" for gn in gs_gi})
    print("Gi^Gq (red)", " ".join(sorted(list(gi_gq), key=lambda gn: Segment.generic_number_of(gn).index)))
    colormap.update({gn: "red" for gn in gi_gq})
    print("Gq^Gs (green)", " ".join(sorted(list(gq_gs), key=lambda gn: Segment.generic_number_of(gn).index)))
    colormap.update({gn: "green" for gn in gq_gs})

    # Fig. 2b
    # Look up these residues in ADRB2
    residues = []
    gen_nums = []
    adrb2 = gpcrdb.GPCRdbEntry('adrb2_human', 'P07550', 'Class A (Rhodopsin)', force=False)
    with open(adrb2.alignment_path) as f:
        for l in f:
            if l.startswith('#'):
                continue
            cols = l.strip().split(',')
            generic_number = cols[2]
            if generic_number == 'None':
                continue
            structure_based_number = generic_number.split('.')[0] + 'x' + generic_number.split('x')[-1]
            if structure_based_number in roi:
                residues.append(cols[0])
                gen_nums.append(generic_number)

    # Write PyMOL commands
    commands = ["fetch 3sn6", "hide everything", "show cartoon, (chain R and resi 1-400) or chain A",
                "color gray90, chain R and elem C", "color lightorange, chain A and elem C", "bg_color white"]

    for r, gen_num in zip(residues, gen_nums):
        commands.append("select {}_{}, chain R and resi {}".format(r, gen_num, r[1:]))
        commands.append("show spheres, {} and name CA".format(r))
        gn = gen_num.split('.')[0] + 'x' + gen_num.split('x')[-1]
        commands.append("color {}, {} and name CA".format(colormap[gn], r))

    commands.append("""set_view (\
     0.945073068,   -0.057450056,    0.321734518,\
    -0.309198350,    0.161730975,    0.937132597,\
    -0.105872914,   -0.985147774,    0.135085553,\
    -0.000358216,    0.000502050, -285.955871582,\
    23.052556992,   10.656435966,   18.293930054,\
   239.106460571,  332.808990479,  -20.000000000 )""")
    commands.append("scene side, store")
    commands.append("hide cartoon, chain A")
    commands.append("""set_view (\
     0.981791675,    0.187008470,    0.032943685,\
    -0.183039993,    0.978211105,   -0.097853623,\
    -0.050524112,    0.090041943,    0.994643092,\
    -0.000335140,   -0.000247810, -157.019104004,\
    15.821222305,    6.886687756,   16.715213776,\
   137.008300781,  176.940826416,  -20.000000000 )""")
    commands.append("scene cavity, store")

    with open("./figures/3b_pymol_commands.pml", 'w') as f:
        f.write('\n'.join(commands))

    # Fig. 2c
    freqs_by_coupling = {primary: [] for primary in GproteinCoupling}
    number_of_receptors = {primary: 0 for primary in GproteinCoupling}
    for receptor in gpcrdb.get_filtered_receptor_list():
        if receptor.receptor_class != 'Class A (Rhodopsin)':
            continue

        number_of_receptors[receptor.primary_coupling] += 1

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
                structure_based_number = anno.generic_number.split('.')[0] + 'x' + latter

                if structure_based_number in roi:
                    freqs_by_coupling[receptor.primary_coupling].append(anno.snv.AF)
    
    print(number_of_receptors)
    print({k: len(v) for k, v in freqs_by_coupling.items()})
    print({k: max(v) for k, v in freqs_by_coupling.items()})

    fig, ax = plt.subplots(1, 1, figsize=(4, 4), dpi=300)

    freqs = [freqs_by_coupling[g] for g in GproteinCoupling]
    colors = [g.color for g in GproteinCoupling]
    labels = [g.value for g in GproteinCoupling]
    ax.hist(freqs, bins=20, range=(0, 1), color=colors, stacked=True, label=labels)
    ax.legend()
    ax.set_xlabel("Allele Freq.")
    ax.set_ylabel("Number of Variants")
    ax.set_yscale('log')

    fig.tight_layout()
    fig.savefig("./figures/3c_contacts.pdf")

if __name__ == '__main__':
    # analyze_positions()
    # analyze_arginine_3x50()
    analyze_G_protein_contact_positions()