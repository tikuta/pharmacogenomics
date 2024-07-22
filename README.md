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
    - AlphaMissense: https://doi.org/10.5281/zenodo.8208688
    - Great Ape Genome Project: https://eichlerlab.gs.washington.edu/greatape/data.html
2. Execute `preprare.py`. The script depends on the following sites:
    - GPCRdb
    - UniProt
    - EnsEMBL
3. Execute `analyze.py`.
