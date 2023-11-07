# Pharmacogenomics for GPCRs

## Prerequirements
### Softwares
You may install via APT system in Ubuntu.
- Python 3
- tabix

#### Libraries for python
- Requests
- Biopython
- Matplotlib
- Pandas

## How to run
1. Place VCF, TBI, and GZ files in `data` directory.
    - 54KJPN: https://jmorp.megabank.tohoku.ac.jp/downloads/tommo-54kjpn-20230626-af_snvindelall
    - 1KGP: https://www.ebi.ac.uk/eva/?eva-study=PRJEB30460
    - AlphaMissense: https://console.cloud.google.com/storage/browser/dm_alphamissense
2. Execute `run.py`. The script depends on the following sites:
    - GPCRdb
    - UniProt
    - EnsEMBL
3. Execute `analyze.py`.
