#!/usr/bin/env python3

import requests
import json

def get_entry(accession, fpath):
    uri = "https://rest.uniprot.org/uniprotkb/{}.json".format(accession)
    r = requests.get(uri)

    if not r.ok:
        raise Exception
    
    j = r.json()
    
    with open(fpath, 'w') as f:
        json.dump(j, f)

    return j

def uniprot2ensembl(entry):
    if entry['uniProtkbId'] == 'GP179_HUMAN':
        # In UniProt, GPR179 is united with ENSG00000276469 on CHR_HSCHR17_7_CTG4.
        # However, GPR179 is mapped on Chromosome 17 in the primary assembly.
        return 'ENSG00000277399'

    for ref in entry['uniProtKBCrossReferences']:
        if ref['database'] == 'Ensembl':
            for prop in ref['properties']:
                if prop['key'] == 'GeneId':
                    gene_id = prop['value'].split('.')[0]
                    return gene_id
    raise Exception("No Gene ID found")