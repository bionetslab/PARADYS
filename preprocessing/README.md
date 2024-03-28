# Preprocessing workflow

We here present our exemplary workflow for generating input files for PARADYS from TCGA PRAD mutation data and SSN networks.



## Necessary files

Download the files `PRAD_nets_SSN.fea`, `TCGA-PRAD_CNVS.tsv`, and `TCGA-PRAD.muse_snv.tsv` from this link. They contain the personalized dysregulation edges from the PRAD cohort, and the CNV and SNV raw mutation data from TCGA PRAD, respectively.



## Producing input files

Run the `preprocessing.py` script with updated input file and output file paths. It creates the necessary `networks.csv` and `mutations.csv` files for PARADYS.



## Already preprocessed PRAD input files

We also offer the resulting preprocessed input files `prad_networks.csv` and `prad_mutations.csv` files for the TCGA PRAD cohort under the same link. You can use those to try out PARADYS on real-life input data. 