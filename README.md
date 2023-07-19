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

## How to run
1. Place VCF files and TBI files in `data` directory.
2. Execute `run.py`. The script depends on the following sites:
    - GPCRdb
    - UniProt
    - EnsEMBL
3. Execute `analyze.py`.