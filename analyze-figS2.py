#!/usr/bin/env python3
import gpcrdb
import matplotlib
matplotlib.use('Agg')
matplotlib.rc('pdf', fonttype=42)
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = "Arial"
import ensembl
from utils import VariationType, Segment
import json

def analyze_nonterminal_regions(filename):
    nonterminal = []
    for receptor in gpcrdb.get_filtered_receptor_list():
        # Generic number system is consistent only within the same class.
        if receptor.receptor_class != 'Class A (Rhodopsin)':
            continue

        with open(receptor.ensembl_path) as f:
            display_name = json.load(f)['display_name']
        with open(receptor.alignment_path) as f:
            seq = ''.join([l.split(',')[0][0] for l in f if not l.startswith('#')])
        with open(receptor.japan_cds_csv_path) as f:
            for l in f:
                try:
                    anno = ensembl.Annotation.from_csv_line(l)  
                except ensembl.BlankLineError:
                    continue
                
                if anno.var_type != VariationType.MISSENSE:
                    continue
                if anno.snv.AF <= 0.5:
                    continue
                if anno.segment == Segment.Nterm or anno.segment == Segment.Cterm:
                    continue

                assert(seq[anno.residue_number - 1] == anno.ref_aa)
                nonterminal.append({"display_name": display_name, "annotation": anno})
    # A.BxC, in many case, B = C
    generic_numbers = set(d['annotation'].generic_number for d in nonterminal if d['annotation'].generic_number)
    # AxC
    structure_based_numbers = set(n.split('.')[0] + 'x' + n.split('x')[-1] for n in generic_numbers)
    
    # Look up these residues in ADRB2
    residues = []
    gen_nums = []
    colors = []
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
            if structure_based_number in structure_based_numbers:
                residues.append(cols[0])
                gen_nums.append(generic_number)
                rgba = Segment.value_of(cols[1]).color
                colors.append("0x{:02x}{:02x}{:02x}".format(int(rgba[0] * 255), int(rgba[1] * 255), int(rgba[2] * 255)))
                print("{} ({})".format(cols[0], structure_based_number))

    # Write PyMOL commands
    commands = ["load https://opm-assets.storage.googleapis.com/pdb/3sn6.pdb", "remove chain A+B+G+N", "hide everything"]
    commands += ["show spheres, resn DUM", "color gray50, resn DUM", "set sphere_scale, 0.1, resn DUM"]
    commands += ["show cartoon, resi 1-400", "color gray90, elem C", "bg_color white"]
    commands.append("""set_view (\
     0.023613820,   -0.006613106,    0.999683261,\
     0.999707758,    0.004089679,   -0.023586458,\
    -0.003932271,    0.999952376,    0.006706934,\
     0.000007592,    0.000207097, -329.484893799,\
    -4.564888000,   -0.986666739,   -5.704017639,\
   249.196655273,  409.797607422,  -20.000000000 )""")
    for r, c, gn in zip(residues, colors, gen_nums):
        commands.append("select {}_{}, resi {}".format(r, gn, r[1:]))
        commands.append("show spheres, {} and name CA".format(r))
        commands.append("color {}, {} and name CA".format(c, r))
    with open(filename, 'w') as f:
        f.write('\n'.join(commands))

def visualize_ptgdr2_v204a(filename):
    commands = ["load https://opm-assets.storage.googleapis.com/pdb/6d26.pdb", "remove resi 1238-1362", "hide everything"]
    commands += ["show spheres, resn DUM", "color gray50, resn DUM", "set sphere_scale, 0.1, resn DUM"]
    commands += ["show cartoon, polymer", "color gray90, elem C", "bg_color white"]
    commands += ["show sticks, resi 204 and not name C+N+O", "color yelloworange, resi 204 and elem C"]
    commands.append("""set_view (\
     0.999601483,    0.023199057,    0.014749675,\
     0.014575680,    0.007802706,   -0.999853611,\
    -0.023312025,    0.999673069,    0.007462042,\
    -0.000005939,    0.000136152, -242.724426270,\
     0.188465744,   -0.155589461,   -2.518058062,\
   162.436721802,  323.037689209,  -20.000000000 )""")
    with open(filename, 'w') as f:
        f.write('\n'.join(commands))

if __name__ == '__main__':
    analyze_nonterminal_regions("./figures/2f_pymol_commands.pml")
    visualize_ptgdr2_v204a("./figures/S2c_pymol_commands.pml")