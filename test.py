import json
import glob

chromosomes = set()
for p in glob.glob("apes/pan_troglodytes/*/*.json"):
    if len(p) > 1:
        print(p)