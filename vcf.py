#!/usr/bin/env python3

import subprocess
import os
from utils import Region
from typing import List

def filter_vcf(vcfs, regions: List[Region], fpath, force=False):
    if os.path.exists(fpath) and force is False:
        return

    with open(fpath, 'w') as f:
        for vcf in vcfs:
            args = ["tabix", vcf] # Run external `tabix` command.
            for region in regions:
                args.append("chr{}:{}-{}".format(region.chromosome, region.start, region.end))

            cp = subprocess.run(args, check=True, text=True, capture_output=True)
            if cp.stdout:
                f.write(cp.stdout)