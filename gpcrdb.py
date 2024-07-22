#!/usr/bin/env python3

import requests
import json
import os
from typing import List, Dict
from config import *
import pandas as pd
from utils import Segment, GproteinCoupling

class GPCRdbEntry:
    def __init__(self, entry_name: str, accession: str, receptor_class: str, primary_coupling=None, force=False) -> None:
        self.entry_name = entry_name
        self.accession = accession
        self.receptor_class = receptor_class
        self.primary_coupling = primary_coupling

        self.dirpath = os.path.join(self.receptor_class, self.entry_name)
        os.makedirs(self.dirpath, exist_ok=True)

        self.generic_number_path = os.path.join(self.dirpath, "gpcrdb.json")
        self.uniprot_path = os.path.join(self.dirpath, "uniprot.json")
        self.alignment_path = os.path.join(self.dirpath, "alignment.csv")
        self.ensembl_path = os.path.join(self.dirpath, "ensembl.json")
        self.sequence_path = os.path.join(self.dirpath, "sequence.json")
        self.cds_path = os.path.join(self.dirpath, "cds.json")
        self.japan_gene_vcf_path = os.path.join(self.dirpath, VCF_JPN_GENE_FILENAME)
        self.japan_cds_vcf_path = os.path.join(self.dirpath, VCF_JPN_CDS_FILENAME)
        self.japan_cds_csv_path = os.path.join(self.dirpath, CSV_JPN_CDS_FILENAME)
        self.global_gene_vcf_path = os.path.join(self.dirpath, VCF_GLOBAL_GENE_FILENAME)
        self.global_cds_vcf_path = os.path.join(self.dirpath, VCF_GLOBAL_CDS_FILENAME)
        self.global_cds_csv_path = os.path.join(self.dirpath, CSV_GLOBAL_CDS_FILENAME)

        self._get_generic_numbers(force=force)
        with open(self.generic_number_path) as f:
            generic_numbers = json.load(f)
            self.protein_seq = ''
            self.generic_numbers = []
            self.segments = []

            for res in generic_numbers:
                self.protein_seq += res['amino_acid']
                self.generic_numbers.append(res['display_generic_number'])
                self.segments.append(Segment.value_of(res['protein_segment']))

        self._get_uniprot_entry(force=force)
        with open(self.uniprot_path) as f:
            j = json.load(f)

            self.ensembl_id_candidates = []
            for ref in j['uniProtKBCrossReferences']:
                if ref['database'] == 'Ensembl':
                    for prop in ref['properties']:
                        if prop['key'] == 'GeneId':
                            gene_id = prop['value'].split('.')[0]
                            self.ensembl_id_candidates.append(gene_id)
        if self.accession == 'Q6PRD1': # gp179_human
            # Overwrite associated ID ENSG00000276469 with manual ID mapping
            self.ensembl_id_candidates = ['ENSG00000277399']
        if len(self.ensembl_id_candidates) == 0:
            # Manual mapping for unassociated entries
            if self.accession == 'P59539': # t2r45_human
                self.ensembl_id_candidates = ['ENSG00000261936']
            elif self.accession == 'Q7Z7M1': # agrd2_human
                self.ensembl_id_candidates = ['ENSG00000180264']
            else:
                raise Exception("No Gene ID found for {}".format(self.entry_name))
        self.ensembl_id_candidates = set(self.ensembl_id_candidates)
        
    def __str__(self) -> str:
        return "{} > {} ({})".format(self.receptor_class, self.entry_name, self.accession)
    
    def _get_generic_numbers(self, force=False):
        if not os.path.exists(self.generic_number_path) or force is True:
            uri = "https://gpcrdb.org/services/residues/extended/{}/".format(self.entry_name)

            r = requests.get(uri)
            if not r.ok:
                raise Exception
            
            j = r.json()
            
            with open(self.generic_number_path, 'w') as f:
                json.dump(j, f, indent=2)

    def _get_uniprot_entry(self, force=False):
        if not os.path.exists(self.uniprot_path) or force is True:
            uri = "https://rest.uniprot.org/uniprotkb/{}.json".format(self.accession)
            r = requests.get(uri)

            if not r.ok:
                raise Exception
            
            j = r.json()
            
            with open(self.uniprot_path, 'w') as f:
                json.dump(j, f, indent=2)

def get_filtered_receptor_list(force=False):
    primary_coupling = _get_primary_coupling(force=force)

    receptors = _get_receptor_list(force=force)
    for r in receptors:
        species = r['species']
        if species != "Homo sapiens":
            continue

        entry_name = r['entry_name']
        if entry_name in BLOCKLIST:
            continue
        
        accession = r['accession']
        receptor_class = r['receptor_class']

        coupling = primary_coupling.get(entry_name.replace("_human", "").upper(), GproteinCoupling.Unknown)

        try:
            yield GPCRdbEntry(entry_name, accession, receptor_class, primary_coupling=coupling, force=force)
        except Exception as e:
            print("Skipped {}: {}".format(entry_name, str(e)))

def _get_receptor_list(force=False) -> List[Dict]:
    p = os.path.join("data", "receptors.json")
    if not os.path.exists(p) or force is True:
        uri = "https://gpcrdb.org/services/receptorlist/"
        r = requests.get(uri)

        if not r.ok:
            raise Exception
        
        j = r.json()
        
        with open(p, 'w') as f:
            json.dump(j, f, indent=2)

    with open(p) as f:
        return json.load(f)

def _get_primary_coupling(force=False):
    html = os.path.join("data", "couplings.html")
    csv = os.path.join("data", "couplings.csv")

    if not os.path.exists(html) or force is True:
        uri = "https://gproteindb.org/signprot/couplings#"
        r = requests.get(uri)

        if not r.ok:
            raise

        with open(html, 'w', encoding='utf-8') as f:
            f.write(r.text)

    dfs = pd.read_html(html, displayed_only=False, attrs={'id': 'familiestabletab'})
    assert(len(dfs) == 1)
    df = dfs[0].drop(0).rename(columns=lambda x: x.replace('  ×', ''))
    df = df[['Receptor', 'Selectivity', 'Family rank orders']][df[('Source', 'Lab')] == 'GproteinDb'].droplevel(0, axis=1)
    df.to_csv(csv)

    df_A = df[df['Cl'] == 'A']

    primary_coupling = {}
    primary_coupling.update({name: GproteinCoupling.Gs for name in set(df_A['Uniprot'][df_A['Primaryfamily'] == "Gs"])})
    primary_coupling.update({name: GproteinCoupling.Gio for name in set(df_A['Uniprot'][df_A['Primaryfamily'] == "Gi/o"])})
    primary_coupling.update({name: GproteinCoupling.Gq11 for name in set(df_A['Uniprot'][df_A['Primaryfamily'] == "Gq/11"])})
    primary_coupling.update({name: GproteinCoupling.G1213 for name in set(df_A['Uniprot'][df_A['Primaryfamily'] == "G12/13"])})

    return primary_coupling