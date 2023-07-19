#!/usr/bin/env python3

import requests
import os
import json

def main():
    save_dir = './uniprotkb'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    with open("../GtoP/GPCRTargets.csv") as f:
        for l in f.readlines()[1:]:
            if len(l) == 0:
                continue
            cols = l.split('","')
            gene = cols[10] # HGNC symbol
            if len(gene) == 0: # some genes are not present in human (only in rodents)
                continue

            swissprot_id = cols[15]
            if len(swissprot_id) == 0: # some genes are pseudogenes
                continue
            print(gene, '->', swissprot_id)
            uri = "https://rest.uniprot.org/uniprotkb/{}.json".format(swissprot_id)
            r = requests.get(uri, headers={"Content-Type": "application/json"})

            if not r.ok:
                raise Exception

            with open(os.path.join(save_dir, gene + '.json'), 'w') as w:
                json.dump(r.json(), w)


if __name__ == '__main__':
    main()