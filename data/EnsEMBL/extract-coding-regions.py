#!/usr/bin/env python3

import os
import json
import glob

with open("coding-regions.txt", 'w') as f:
    for p in glob.glob("./lookup/*.json"):
        gene = os.path.splitext(os.path.basename(p))[0]

        j = json.load(open(p, 'r'))

        biotype = j['biotype']
        if biotype != "protein_coding":
            continue

        chr_num, strand = j['seq_region_name'], j['strand']
        
        transcripts = j['Transcript']
        for t in transcripts:
            is_canonical = t['is_canonical'] # 0 or 1
            if int(is_canonical) == 0: # we consider canonical isoforms only
                continue

            assert('Translation' in t) # "canonical" assumes proteins

            ensp_id = t['Translation']['id']

            # We only have positions where coding regions start and end
            #                    and where exons          start and end.
            # So we must assembly exact coding regions (indicated as '*' in the following examples).
            # Note that start (s) is smaller than end (e), reagrdless of strand (plus or minus).
            # 
            #                   cds_s                                cds_e
            #                     |                                    |
            # ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
            #                     *****       *********       **********
            #  |---|     |------------|       |-------|       |--------------|        |---| 
            #  exon1   exon2_s      exon2_e exon3_s exon3_e exon4_s        exon4_e    exon5
            #
            # Or,
            #                   cds_s                                 cds_e
            #                     |                                     |
            # ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
            #                     ***************************************
            #             |------------------------------------------------------------|
            #           exon6_s                                                      exon6_e

            cds_s = t['Translation']['start']
            cds_e = t['Translation']['end']
            assert(cds_s < cds_e)

            coding_regions = []
            for exon in t['Exon']:
                exon_s = int(exon['start'])
                exon_e = int(exon['end'])
                assert(exon_s < exon_e)
                
                if exon_s <= cds_s <= exon_e <= cds_e:   # Exon 2 in the above example
                    coding_regions.append([cds_s, exon_e])
                elif cds_s <= exon_s <= exon_e <= cds_e: # Exon 3
                    coding_regions.append([exon_s, exon_e])
                elif cds_s <= exon_s <= cds_e <= exon_e: # Exon 4
                    coding_regions.append([exon_s, cds_e])
                elif exon_s <= cds_s <= cds_e <= exon_e: # Exon 6
                    coding_regions.append([cds_s, cds_e])
                else:                                    # Exons 1 and 5 contain no coding region
                    continue

            cols = [gene, ensp_id, "chr" + chr_num, str(strand)]
            cols += ["{}-{}".format(r[0], r[1]) for r in coding_regions]
            f.write('\t'.join(cols) + '\n')
