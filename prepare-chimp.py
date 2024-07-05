import gpcrdb
import requests
import os
import time
import json

class EnsemblHomology:
    def __init__(self, gpcrdb_entry: gpcrdb.GPCRdbEntry, force=False) -> None:
        self.gpcrdb_entry = gpcrdb_entry
        self.homolog_path = os.path.join(gpcrdb_entry.dirpath, "homology.json")
        self._get_homologs(force=force)

        with open(self.homolog_path) as f:
            j = json.load(f)
            if not (1 == len(j['data']) == len(j['data'][0]['homologies'])):
                print("Error on", gpcrdb_entry)
                return
            self.protein_id = j['data'][0]['homologies'][0]['target']['protein_id']
            self.gene_id = j['data'][0]['homologies'][0]['target']['id']

    def _get_homologs(self, force=False, max_try=5):
        if os.path.exists(self.homolog_path) and force is False:
            print("Skipping", self.gpcrdb_entry.entry_name)
            return

        to_dump = None
        for gene_id in self.gpcrdb_entry.ensembl_id_candidates:
            # pan_troglodytes = chimpanzee
            uri = 'https://rest.ensembl.org/homology/id/human/{}?target_species=pan_troglodytes'.format(gene_id)

            r = requests.get(uri, headers={"Content-Type": "application/json"})

            if not r.ok:
                try_count = 1
                retry_waits = [None, 10, 30, 60, 120, 300]
                while try_count <= max_try:
                    print("Request to EnsEMBL lookup API failed (try = {}). Waiting for {} seconds...".format(try_count, retry_waits[try_count]))
                    time.sleep(retry_waits[try_count])
                    r = requests.get(uri, headers={"Content-Type": "application/json"})
                    if r.ok:
                        break
                    try_count += 1
                else:
                    raise Exception
            
            j = r.json()

            if to_dump is None:
                to_dump = j
            with open(self.homolog_path, 'w') as f:
                json.dump(to_dump, f, indent=2)
        

class EnsemblChimpanzeeGene:
    def __init__(self, gene_id: str, gpcrdb_entry: gpcrdb.GPCRdbEntry, force=False) -> None:
        self.gene_id = gene_id
        self.gpcrdb_entry = gpcrdb_entry
        self.gene_path = os.path.join(self.gpcrdb_entry.dirpath, "chimpanzee.json")
        self._get_gene(force=force)

    def _get_gene(self, force=False, max_try=5):
        if os.path.exists(self.gene_path) and force is False:
            print("Skipping", self.gpcrdb_entry.entry_name)
            return

        to_dump = None
        for gene_id in self.gpcrdb_entry.ensembl_id_candidates:
            uri = 'https://rest.ensembl.org/lookup/id/{}?expand=1'.format(gene_id)

            r = requests.get(uri, headers={"Content-Type": "application/json"})

            if not r.ok:
                try_count = 1
                retry_waits = [None, 10, 30, 60, 120, 300]
                while try_count <= max_try:
                    print("Request to EnsEMBL lookup API failed (try = {}). Waiting for {} seconds...".format(try_count, retry_waits[try_count]))
                    time.sleep(retry_waits[try_count])
                    r = requests.get(uri, headers={"Content-Type": "application/json"})
                    if r.ok:
                        break
                    try_count += 1
                else:
                    raise Exception
            
            j = r.json()

            if to_dump is None:
                to_dump = j
            with open(self.gene_path, 'w') as f:
                json.dump(to_dump, f, indent=2)

def main():
    for gpcrdb_entry in gpcrdb.get_filtered_receptor_list():
        homology = EnsemblHomology(gpcrdb_entry)
        homolog_gene = EnsemblChimpanzeeGene(homology.gene_id, gpcrdb_entry)



if __name__ == '__main__':
    main()