#!usr/bin/env python3

import numpy as np
import glob
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = "Arial"

def is_snv(ref, alt):
    if len(ref) == len(alt) == 1:
        return True
    if len(ref) > 1:
        return False
    if len(alt) > 1:
        if ',' in alt:
            bases = alt.split(',')
            if 1 in [len(b) for b in bases]:
                return True
        else:
            return False

def curate_SNV_sites(vcf):
    all = set()
    with open(vcf) as f:
        snv = set()
        for l in f.readlines():
            if len(l) > 0:
                cols = l.split('\t')
                chr, pos, ref, alt = cols[0], cols[1], cols[3], cols[4]
                all.add(chr + ":" + pos)
                if is_snv(ref, alt):
                    snv.add(chr + ":" + pos)
        return snv, all

def analyze():
    all_sites = set()

    sites_coding = set()
    for vcf in glob.glob("../data/38KJPN/coding-regions/*.vcf"):
        snv, _ = curate_SNV_sites(vcf)
        sites_coding |= snv

    sites_gene = set()
    for vcf in glob.glob("../data/38KJPN/gene-regions/*.vcf"):
        snv, all = curate_SNV_sites(vcf)
        sites_gene |= snv
        all_sites |= all

    assert(sites_coding.issubset(sites_gene))
    assert(sites_gene.issubset(all_sites))

    return len(sites_coding), len(sites_gene), len(all_sites)

def plot_SNV_calls_on_GPCR_coding_regions(num_coding, num_gene):
    fig, ax = plt.subplots(1, 1, figsize=(4, 2), dpi=300)
    ax.barh(0, num_coding, color='tab:orange')
    ax.barh(0, num_gene - num_coding, left=num_coding, color='tab:gray', alpha=0.6)
    ax.set_xlim(0, num_gene)
    non_coding_text = "Non-coding region\n{:,} calls ({:.1f}%)".format(num_gene - num_coding, (num_gene - num_coding) / num_gene * 100)
    ax.text((num_gene - num_coding) / 2, 0, non_coding_text, ha='center', va='center')
    coding_text = "Coding region\n{:,} calls ({:.1f}%)".format(num_coding, num_coding / num_gene * 100)
    ax.text(num_coding / 2, -0.45, coding_text, ha='center', va='top')
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig("SNV_calls_on_GPCR_coding_regions.pdf")

def main():
    num_coding, num_gene, num_sites = analyze() # 35974, 1248615, 1397696
    print(num_sites)

    plot_SNV_calls_on_GPCR_coding_regions(num_coding, num_gene)

if __name__ == '__main__':
    main()
