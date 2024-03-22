# Overview

PAtient-specific RAnking of genes driving DYSregulations

## Installation
To use PARADYS, you need to have Python installed on your system.
Install conda environment as follows
```
conda env create -f environment.yml
```

## Usage
You can use PARADYS from the command line with the following arguments:
```
python paradys.py patient networks mutations kappa rho d scores
```
Where:
* patient: patient ID (str)
* networks: path to the dysregulation networks file (str)
* mutations: path to the mutation matrix file (str)
* kappa: parameter kappa (int) 
* rho: parameter rho (optional, default is 0).
* d: dumping factor (optional, default is 0.85).
* scores: specify if you want to return scores (optional, default is False)

Example usage:
```
python paradys.py 'TCGA-AA-0000' network_file.csv mutation_matrix.csv 3 0 0.85 True
```

## Output
The patient data will be processed using the specified parameters and provide the predictions (and dysregulation scores) as output.


## Evaluating PARADYS
For a large-scale empirical evaluation of PARADYS, please follow the instructions given here: 
https://github.com/bionetslab/PARADYS/tree/main/validation

