# PARADYS

In this repository we present our tool PARADYS for **PA**tient-specific **RA**nking of genes driving **DYS**regulations in cancer. Below you can find information on installation, usage and exemplary calls.



## Installation

To use PARADYS you need to have Python installed on your system. Install the necessary conda environment by running

```
conda env create -f environment.yml
```

and then activate the newly created environment.



## Usage

You can use PARADYS from the command line with the following arguments:
```
python paradys.py --patients <PATIENT_IDS> --networks <NETWORK_FILE> --mutations <MUTATION_FILE> --outdir <OUTPUT_DIR> --kappa <KAPPA> --d <D> [--scores] [--directed]
```
The input arguments are given by:
* `PATIENT_IDS`: List of space-separated patient IDs of patients that are supposed to be analyzed. The patients IDs must match with those used in your input files.
* `NETWORK_FILE`: Path to the input dysregulation networks file.
* `MUTATION_FILE`: Path to the input mutation matrix file.
* `OUTPUT_DIR`: Path to output directory for storing results.
* `KAPPA`: Value of parameter kappa (int). Sets size of considered shortest-path neighborhood for putative driver gene identification. Larger values correspond to potentially more detected drivers, but also increase runtime.
* `D`: Dumping factor for PageRank algorithm (optional, default is 0.85). Only applies if `--scores` is set.
* `--scores`: Set this flag if you want to compute impact scores for drivers. (optional, default is False)
* `--directed`: Set this flag if your input dysregulation networks are directed. This will be accounted for in the ranking process. Only applies if `--scores` is set. (optional, default is False).



## Input

PARADYS expects two comma-separated `.csv` input files: 

1) The input `NETWORK_FILE` must contain the columns `'patient'`, `'tf'`, and `'gene'`. It has to contain all directed dysregulation edges `(a,b)` from your personalized input networks, with the source node `a` stored in column `'tf'` and target node `b` stored in `'gene'`. The column `'patient'` stores the ID of the patient that this edge belongs to. 
2) The input `MUTATION_FILE` must contain the two columns `'patient'` and `'gene'`. It has to contain the sets of mutated genes for each patient, with the mutated gene in the column `'gene'` and the corresponding patient ID in the column `'patient'`.



## Example Call

We added two toy input data files under `example/networks_mini.csv` and `examples/mutations_mini.csv`. For analyzing patients `P2` and `P3` you can then run PARADYS on this input by calling:

```
python paradys.py --patients P2 P3 --networks examples/networks_mini.csv --mutations examples/mutations_mini.csv --outdir results/ --kappa 2
```

For additionally computing per-patient driver impact scores you can simply run
```
python paradys.py --patients P2 P3 --networks examples/networks_mini.csv --mutations examples/mutations_mini.csv --outdir results/ --kappa 2 --scores
```

If you want to use PARADYS on "real-life" input data, follow the download instructions in `preprocessing/README.md`. Be aware that this input cohort consists of ~470 patients and therefore takes significantly longer to compute. 

## Output

In the output directory, for each analyzed patient we create a file named `<PATIENT_ID>_drivers.csv` storing detected significant driver genes for each analyzed patient. In the column `'driver'` we store the detected significant driver gene, in combination with the associated dysregulation edge in the column `'dysregulation'` and the respective p-value of the chi-squared test in the column `'pvalue'`.

If additionally the `--scores` flag was set, we create a per-patient score file named `<PATIENT_ID>_scores.csv` storing drivers in the column `'driver'` and the corresponding computed impact score in `'score'`.
