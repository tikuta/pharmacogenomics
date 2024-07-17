import gpcrdb
import os
import json
import enum
from utils import GET_json_with_retries

class GreatApe(enum.Enum):
    # human = "homo_sapiens"
    chimpanzee = "pan_troglodytes"
    bonobo = "pan_paniscus"
    gorilla = "gorilla_gorilla"
    sumatran_orangutan = "pongo_abelii"
    # bornean_orangutan = "pongo_pygmaeus" # No entries in EnsEMBL database


class EnsemblHomology:
    def __init__(self, gpcrdb_entry: gpcrdb.GPCRdbEntry, force=False) -> None:
        self.gpcrdb_entry = gpcrdb_entry
        self.dirnames = {ape: os.path.join("apes", ape.value, gpcrdb_entry.entry_name) for ape in GreatApe}
        for d in self.dirnames.values():
            os.makedirs(d, exist_ok=True)
        self._get_homologs(force=force)

        protein_ids = {ape: [] for ape in GreatApe}
        gene_ids = {ape: [] for ape in GreatApe}
        for ape in GreatApe:
            for gene_id in self.gpcrdb_entry.ensembl_id_candidates:
                p = os.path.join(self.dirnames[ape], f"{gene_id}.json")
                with open(p) as f:
                    j = json.load(f)

                    assert(len(j['data']) == 1)

                    for h in j['data'][0]['homologies']:
                        protein_ids[ape].append(h['target']['protein_id'])
                        gene_ids[ape].append(h['target']['id'])
        self.protein_ids = protein_ids
        self.gene_ids = gene_ids

    def _get_homologs(self, force=False):
        for ape in GreatApe:
            for gene_id in self.gpcrdb_entry.ensembl_id_candidates:
                p = os.path.join(self.dirnames[ape], f"{gene_id}.json")
                if os.path.exists(p) and force is False:
                    continue

                uri = f'https://rest.ensembl.org/homology/id/human/{gene_id}?target_species={ape.value}'
                j = GET_json_with_retries(uri)
            
                with open(p, 'w') as f:
                    json.dump(j, f, indent=2)

def main():
    num_homologs = {ape: 0 for ape in GreatApe}
    for gpcrdb_entry in gpcrdb.get_filtered_receptor_list():
        homology = EnsemblHomology(gpcrdb_entry)
        for ape in GreatApe:
            num_homologs[ape] += len(homology.protein_ids[ape])

    print(num_homologs)

if __name__ == '__main__':
    main()