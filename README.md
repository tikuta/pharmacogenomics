# Pharmacogenomics for GPCRs
## Prerequirements
### Softwares
- Python 3.8.10
    - Requests
    - Biopython
    - Matplotlib
    - Pandas
    - (optional) CrossMap 0.7.0
- tabix 1.10.2
- (optional) samtools 1.10.2
- (optional) bcftools 1.10.2

## How to run
1. Execute `wget.sh` to download files. 
    - 54KJPN (`./data/54KJPN`)
        - https://jmorp.megabank.tohoku.ac.jp/downloads/tommo-54kjpn-20230626-af_snvindelall
    - 1KGP (`./data/1KGP`)
        - https://www.ebi.ac.uk/eva/?eva-study=PRJEB30460
    - AlphaMissense (`./data`)
        - https://doi.org/10.5281/zenodo.8208688
    - (optional) Great Ape Genome Project (`./data/GAGP`)
        - https://eichlerlab.gs.washington.edu/greatape/data.html
        - https://hgdownload.soe.ucsc.edu/downloads.html#chimp
        - GAGP dataset needs to be remapped on panTro5.
            - `gzip -d -f -k panTro5.fa.gz`
            - `samtools faidx -o panTro5.fa.fai panTro5.fa`
            - `CrossMap vcf --compress panTro2ToPanTro5.over.chain.gz Pan_troglodytes.vcf.gz panTro5.fa crossmap.vcf.gz`
            - `bcftools sort -Oz -o crossmap_sorted.vcf.gz crossmap.vcf.gz`
            - `tabix -p vcf crossmap_sorted.vcf.gz`
2. Execute `preprare.py`. The script depends on the following sites:
    - GPCRdb
        - https://gpcrdb.org/
    - UniProt
        - https://www.uniprot.org/
    - EnsEMBL
        - https://useast.ensembl.org/index.html
3. Execute `analyze-*.py` and `analyze-*.ipynb`.
