#!/usr/bin/env python3
import gpcrdb
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = "Arial"
import ensembl
from utils import VariationType, Segment, GproteinCoupling
import json

def analyze_positions():
    assigned = {}
    found = {}
    for receptor in gpcrdb.get_filtered_receptor_list():
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

    family_stats = {}
    aa_stats = {}
    codon_stats = {}
    for receptor in gpcrdb.get_filtered_receptor_list():
        family_stats[receptor.receptor_class] = family_stats.get(receptor.receptor_class, 0) + 1

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
    
    fig, ax = plt.subplots(1, 1, figsize=(8, 2), dpi=300)

    left = 0
    for family in sorted(family_stats.keys()):
        delta = family_stats[family]
        color = 'tab:orange' if family == 'Class A (Rhodopsin)' else 'lightgray'
        ax.barh(2, delta, left=left, color=color, linewidth=0.5, edgecolor='k')
        if family == 'Class A (Rhodopsin)':
            ax.text(left + delta / 2, 2, "Family A\n({})".format(delta), ha='center', va='center', size=6)
        elif family.startswith("Class"):
            family_name = family.split(" ")[1]
            ax.text(left + delta / 2, 2, "{}\n({})".format(family_name, delta), ha='center', va='center', size=6)
        left += delta
    total_receptors = left

    left = 0
    for aa in sorted(aa_stats.keys(), key=lambda aa:aa_stats[aa], reverse=True):
        delta = aa_stats[aa]
        color = 'tab:orange' if aa == 'R' else 'lightgray'
        ax.barh(1, delta, left=left, color=color, linewidth=0.5, edgecolor='k')
        if aa == 'R':
            ax.text(left + delta / 2, 1, "Arg\n({})".format(delta), ha='center', va='center', size=6)
        left += delta
    
    left = 0
    for codon in sorted(codon_stats.keys(), key=lambda codon: ('CG' not in codon, -codon_stats[codon])):
        delta = codon_stats[codon]
        color = 'tab:orange' if 'CG' in codon else 'lightgray'
        ax.barh(0, delta, left=left, color=color, linewidth=0.5, edgecolor='k')
        ax.text(left + delta / 2, 0, "{}\n({})".format(codon, delta), ha='center', va='center', size=6)
        left += delta

    ax.set_xlim(0, total_receptors)
    ax.set_xticks([total_receptors], minor=True)
    ax.set_xticklabels([total_receptors], minor=True)
    ax.set_xlabel("Number of GPCRs")

    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(['3x50 Codon', '3x50 AA', "GPCRs"])
    
    fig.tight_layout()
    fig.savefig("./figures/S3a_arginine_3x50.pdf")

def analyze_G_protein_contact_positions():
    # See the following reference for detail
    # https://doi.org/10.1038/s41467-022-34055-5
    common_residues = {"3x50", "3x53", "3x54", "34x50", "34x51", "34x55", "5x65", "5x68", "6x32", "6x33", "6x36", "6x37", "7x56", "8x47"}
    gs_residues = common_residues | {"3x54", "3x55", "34x51", "34x54", "34x55", "5x64", "5x68", "5x69", "5x71", "5x72", "5x74", "5x75", "5x77", "8x48"}
    gi_residues = common_residues | {"12x49", "2x40", "3x50", "3x53", "34x52", "34x55", "5x71", "6x32", "7x56", "8x47", "8x49"}
    gq_residues = common_residues | {"2x37", "2x39", "2x40", "3x49", "34x51", "34x53", "34x55", "34x56", "34x57", "4x38", "4x39", "6x30", "6x33", "8x48", "8x49"}

    f = lambda gn: (Segment.generic_number_of(gn).index, int(gn.split('x')[-1]))

    roi = frozenset(gs_residues | gi_residues | gq_residues)
    gs_only = frozenset(gs_residues - gi_residues - gq_residues)
    print("Gs only", " ".join(sorted(list(gs_only), key=f)))
    gi_only = frozenset(gi_residues - gs_residues - gq_residues)
    print("Gi only", " ".join(sorted(list(gi_only), key=f)))
    gq_only = frozenset(gq_residues - gs_residues - gi_residues)
    print("Gq only", " ".join(sorted(list(gq_only), key=f)))
    gs_gi = frozenset(gs_residues & gi_residues - gq_residues)
    print("Gs^Gi", " ".join(sorted(list(gs_gi), key=f)))
    gi_gq = frozenset(gi_residues & gq_residues - gs_residues)
    print("Gi^Gq", " ".join(sorted(list(gi_gq), key=f)))
    gq_gs = frozenset(gq_residues & gs_residues - gi_residues)
    print("Gq^Gs", " ".join(sorted(list(gq_gs), key=f)))
    gs_gi_gq = frozenset(gs_residues & gi_residues & gq_residues)
    print("Gs^Gi^Gq", " ".join(sorted(list(gs_gi_gq), key=f)))

    # CMYK coloring (C = Gs, M = Gi, Y = Gq, K = common)
    pymol_colormap = {gs_gi_gq: "gray50",
        gs_only: "cyan", gs_gi: "blue",
        gi_only: "magenta", gi_gq: "red",
        gq_only: "yellow", gq_gs: "green"
    }
    matplot_colormap = {gs_gi_gq: "tab:gray",
        gs_only: GproteinCoupling.Gs.color, gs_gi: "blue",
        gi_only: GproteinCoupling.Gio.color, gi_gq: "red",
        gq_only: GproteinCoupling.Gq11.color, gq_gs: "lime"
    }

    # Fig. 3b
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
        color = next((pymol_colormap[residues] for residues in pymol_colormap.keys() if gn in residues), None)
        commands.append("color {}, {} and name CA".format(color, r))

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

    anno_by_coupling = {primary: {gn: [] for gn in roi} for primary in GproteinCoupling}
    number_of_receptors = {primary: 0 for primary in GproteinCoupling}
    for receptor in gpcrdb.get_filtered_receptor_list():
        if receptor.receptor_class != 'Class A (Rhodopsin)':
            continue

        with open(receptor.ensembl_path) as f:
            display_name = json.load(f)['display_name']

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
                    d = {"annotation": anno, "display_name": display_name}
                    anno_by_coupling[receptor.primary_coupling][structure_based_number].append(d)

    # Fig. 3c
    fig, ax = plt.subplots(1, 1, figsize=(4, 2.5), dpi=300)

    data = [sum([[d['annotation'].snv.AF for d in annos] for gn, annos in anno_by_coupling[g].items() if gn in gs_gi_gq], []) for g in GproteinCoupling]
    labels = [g.value for g in GproteinCoupling]
    colors = [g.color for g in GproteinCoupling]
    ax.hist(data, bins=20, range=(0, 1), label=labels, color=colors, stacked=True, edgecolor='k', lw=0.5, orientation='horizontal')
    for g in GproteinCoupling:
        for gn, annos in anno_by_coupling[g].items():
            if gn not in gs_gi_gq:
                continue
            for d in annos:
                anno = d['annotation']
                if anno.snv.AF > 0.2:
                    y = 0.975 if d['display_name'] == 'GPR148' else 0.225
                    text = " {} {}{}$^{{{}}}${}".format(d['display_name'], anno.ref_aa, anno.residue_number, anno.generic_number, anno.alt_aa)
                    ax.text(1, y, text, ha='left', va='center')
    ax.set_xscale('log')
    ax.set_xlabel("Number of Variants")
    ax.set_ylabel("Allele Freq.")
    ax.legend(title="Primary coupling")

    fig.tight_layout()
    fig.savefig("./figures/3c_contacts.pdf")
    plt.close(fig)       

    # Fig. S3c
    fig, axes = plt.subplots(len(GproteinCoupling), 1, figsize=(5, 5), dpi=300)

    for ax, g in zip(axes, GproteinCoupling):
        left = 0
        for residues in matplot_colormap.keys():
            delta = sum([len(anno_by_coupling[g][gn]) for gn in residues])
            ax.barh(0, delta, left=left, color=matplot_colormap[residues])
            left += delta
        ax.set_xlim(0, left)
        ax2 = ax.twiny()
        ax2.set_xticks([left])
        ax2.set_xticklabels([left])
        ax.set_yticks([0])
        ax.set_yticklabels(["{}\n(n = {})".format(g.value, number_of_receptors[g])])
        if g == GproteinCoupling.Gq11:
            ax.set_ylabel("Primary coupling")
    axes[-1].set_xlabel("Number of Variants")

    fig.tight_layout()
    fig.savefig("./figures/S3c_contacts.pdf")

if __name__ == '__main__':
    analyze_positions()
    analyze_arginine_3x50()
    analyze_G_protein_contact_positions()