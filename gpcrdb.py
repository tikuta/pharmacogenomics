#!/usr/bin/env python3

import requests
import json
import os
from typing import List, Dict
from Bio.Align import PairwiseAligner
from misc import BLOCKLIST
import pandas as pd

def get_filtered_receptor_list(fpath, force=False):
    receptors = get_receptor_list(fpath, force=force)
    for receptor in receptors:
        species = receptor['species']
        if species != "Homo sapiens":
            continue
        entry_name = receptor['entry_name']
        if entry_name in BLOCKLIST:
            continue
        yield receptor

def get_receptor_list(fpath, force=False):
    if os.path.exists(fpath) and force is False:
        with open(fpath) as f:
            return json.load(f)

    uri = "https://gpcrdb.org/services/receptorlist/"
    r = requests.get(uri)

    if not r.ok:
        raise Exception
    
    j = r.json()
    
    with open(fpath, 'w') as f:
        json.dump(j, f, indent=2)

    return j

def get_generic_number(entry_name, fpath, force=False):
    if os.path.exists(fpath) and force is False:
        with open(fpath) as f:
            return json.load(f)

    uri = "https://gpcrdb.org/services/residues/extended/{}/".format(entry_name)

    r = requests.get(uri)
    if not r.ok:
        raise Exception
    
    j = r.json()
    
    with open(fpath, 'w') as f:
        json.dump(j, f, indent=2)
    
    return j

def extract_sequence(generic_numbers) -> str:
    amino_acids = [None] * len(generic_numbers)
    for residue in generic_numbers:
        sequence_number = int(residue['sequence_number'])
        amino_acids[sequence_number - 1] = residue['amino_acid']
    return ''.join(amino_acids)

def extract_segments(generic_numbers) -> List[str]:
    ret = [None] * len(generic_numbers)
    for d in generic_numbers:
        sequence_num = d['sequence_number']
        segment = d['protein_segment']
        ret[sequence_num - 1] = segment
    assert(None not in ret)
    return ret

def extract_generic_numbers(generic_numbers, scheme) -> List[str]:
    ret = [None] * len(generic_numbers)
    for d in generic_numbers:
        sequence_num = d['sequence_number']
        if d['display_generic_number'] is None:
            ret[sequence_num - 1] = None
            continue
        for alt in d['alternative_generic_numbers']:
            if alt['scheme'] == scheme:
                ret[sequence_num - 1] = alt['label']
    return ret

def match(seq, generic_numbers, receptor_class, fpath, alignment_for_human=None, force=False):
    if os.path.exists(fpath) and force is False:
        with open(fpath) as f:
            return json.load(f)
    
    class2scheme = {
        "Class A (Rhodopsin)": "BW", 
        "Class B1 (Secretin)": "Wootten", 
        "Class B2 (Adhesion)": "Wootten", 
        "Class C (Glutamate)": "Pin",
        "Class F (Frizzled)": "Wang",
        "Class T (Taste 2)": None,
        "Other GPCRs": None
    }
    scheme = class2scheme[receptor_class]
    segments = extract_segments(generic_numbers)
    generic_nums = extract_generic_numbers(generic_numbers, scheme)
    ref_seq = extract_sequence(generic_numbers)

    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.1
    alignment = aligner.align(ref_seq, seq)[0]

    residues = []
    for target, query in zip(*alignment.aligned):
        for t_res_idx, q_res_idx in zip(range(*target), range(*query)):
            residue = {
                "ensembl_sequence_number": q_res_idx + 1, 
                "ensembl_amino_acid": seq[q_res_idx], 
                "segment": segments[t_res_idx], 
                "generic_number": generic_nums[t_res_idx],
                "gpcrdb_sequence_number": t_res_idx + 1,
                "gpcrdb_amino_acid": ref_seq[t_res_idx]
            }
            residues.append(residue)
    ensembl_seq_nums = [r['ensembl_sequence_number'] for r in residues]
    for i in range(len(seq)):
        res_num = i + 1
        if res_num not in ensembl_seq_nums:
            residue = {
                "ensembl_sequence_number": res_num, 
                "ensembl_amino_acid": seq[i], 
                "segment": None, 
                "generic_number": None,
                "gpcrdb_sequence_number": None,
                "gpcrdb_amino_acid": None
            }
            residues.append(residue)
    residues.sort(key=lambda r: r['ensembl_sequence_number'])

    # Fill missing segment
    for i, r in enumerate(residues):
        if r['segment'] is None:
            n_seg = 'N-term'
            for rn in reversed(residues[:i]):
                if rn['segment'] is not None:
                    n_seg = rn['segment']
                    break

            c_seg = 'C-term'
            for rc in residues[i:]:
                if rc['segment'] is not None:
                    c_seg = rc['segment']
                    break
            
            if n_seg == c_seg:
                new_r = r
                new_r['segment'] = n_seg
                residues[i] = new_r

    ret = {"receptor_class": receptor_class, "scheme": scheme, "residues": residues}

    with open(fpath, 'w') as f:
        json.dump(ret, f, indent=2)

    if alignment_for_human:
        with open(alignment_for_human, 'w') as f:
            f.write(str(alignment))

    return ret

def extract_coupling(fpath, force=False) -> pd.DataFrame:
    if os.path.exists(fpath) and force is False:
        return pd.read_csv(fpath)
    
    html = "couplings.html"
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
    df.to_csv("couplings.csv")
    extracted = df[['Source', 'Receptor', 'Guide to Pharmacology']][df[('Source', 'Group')] == 'Inoue'].droplevel(0, axis=1)
    extracted.to_csv(fpath)

    return pd.read_csv(fpath)

def primary_coupled_receptors() -> Dict:
    df = extract_coupling("couplings_extracted.csv")
    df_A = df[df['Cl'] == 'A']
    df_primary_Gs    = df_A[(df_A['Gs'] == "1'") & (df_A['Gi/o'] != "1'") & (df_A['Gq/11'] != "1'") & (df_A['G12/13'] != "1'")]
    df_primary_Gio   = df_A[(df_A['Gs'] != "1'") & (df_A['Gi/o'] == "1'") & (df_A['Gq/11'] != "1'") & (df_A['G12/13'] != "1'")]
    df_primary_Gq11  = df_A[(df_A['Gs'] != "1'") & (df_A['Gi/o'] != "1'") & (df_A['Gq/11'] == "1'") & (df_A['G12/13'] != "1'")]
    df_primary_G1213 = df_A[(df_A['Gs'] != "1'") & (df_A['Gi/o'] != "1'") & (df_A['Gq/11'] != "1'") & (df_A['G12/13'] == "1'")]
    df_promiscuous   = df_A[~df_A.index.isin(df_primary_Gs.index) &
                      ~df_A.index.isin(df_primary_Gio.index) &
                      ~df_A.index.isin(df_primary_Gq11.index) &
                      ~df_A.index.isin(df_primary_G1213.index)]
    print(df_A.shape, df_primary_Gs.shape, df_primary_Gio.shape, df_primary_Gq11.shape, df_primary_G1213.shape, df_promiscuous.shape)

    d = {
        "Gs":          list(df_primary_Gs['Uniprot']),
        "Gi/o":        list(df_primary_Gio['Uniprot']),
        "Gq/11":       list(df_primary_Gq11['Uniprot']),
        "G12/13":      list(df_primary_G1213['Uniprot']),
        "promiscuous": list(df_promiscuous['Uniprot']),
        }
    return d
    

if __name__ == '__main__':
    pass