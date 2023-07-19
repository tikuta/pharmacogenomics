#!usr/bin/env python3

import numpy as np
import glob
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = "Arial"

def analyze():
    num_missense, num_silent, num_nonsense = 0, 0, 0
    for txt in glob.glob("../data/38KJPN/coding-variations/*.txt"):
        with open(txt) as f:
            print(txt)
            for l in f.readlines():
                if len(l) == 0:
                    continue
                cols = l.split('\t')
                aa_ref = cols[0]
                aa_alts = cols[3].split(',')

                for aa_alt in aa_alts:
                    if aa_alt == '*':
                        num_nonsense += 1
                    elif aa_ref == aa_alt:
                        num_silent += 1
                    else:
                        num_missense += 1
    return num_missense, num_silent, num_nonsense

def plot_SNV_types(num_missense, num_silent, num_nonsense):
    fig, ax = plt.subplots(1, 1, figsize=(4, 2), dpi=300)
    total = num_missense + num_silent + num_nonsense
    ax.barh(0, num_missense, color='tab:orange')
    ax.barh(0, num_silent, left=num_missense, color='tab:gray', alpha=0.6)
    ax.barh(0, num_nonsense, left=num_missense + num_silent, color='tab:gray')
    ax.set_xlim(0, total)

    missense_text = "Missense\n{:,} variations\n({:.1f}%)".format(num_missense, num_missense / total * 100)
    ax.text(num_missense / 2, 0, missense_text, ha='center', va='center')

    silent_text = "Silent\n{:,} variations\n({:.1f}%)".format(num_silent, num_silent / total * 100)
    ax.text(num_silent / 2 + num_missense, 0, silent_text, ha='center', va='center')

    nonsense_text = "Nonsense\n{:,} variations\n({:.1f}%)".format(num_nonsense, num_nonsense / total * 100)
    ax.text(num_nonsense / 2 + num_missense + num_silent, -0.45, nonsense_text, ha='center', va='top')

    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig("SNV_types.pdf")

def main():
    num_missense, num_silent, num_nonsense = analyze() # 23477, 14399, 929
    print(num_missense, num_silent, num_nonsense)
    plot_SNV_types(num_missense, num_silent, num_nonsense)

if __name__ == '__main__':
    main()