import json
import glob

chromosomes = set()
for p in glob.glob("apes/pan_troglodytes/*.json"):
    with open(p) as f:
        j = json.load(f)
        chromosomes.add(j['seq_region_name'])
print(chromosomes)