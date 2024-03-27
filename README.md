# PARADYS

In this repository we present our tool PARADYS for **PA**tient-specific **RA**nking of genes driving **DYS**regulations in cancer. Below you can find information on installation, usage and exemplary calls.



## Installation

To use PARADYS you need to have Python installed on your system. Install the necessary conda environment as follows:

```
conda env create -f environment.yml
```



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
* `KAPPA`: Value of parameter kappa (int). Sets size of considered shortest-path neighborhood for putative driver gene identification. 
* `D`: Dumping factor for PageRank algorithm (optional, default is 0.85).
* `--scores`: Set this flag if you want to compute impact scores for drivers. (optional, default is False)
* `--directed`: Set this flag if your input dysregulation networks are directed. This will be accounted for in the ranking process. Is only used when `--scores` is also set (optional, default is False).



## Example Call

We added two toy input data files under `example/networks_mini.csv` and `examples/mutations_mini.csv`. For analyzing patients `P2` and `P3` you can then run PARADYS on this input by calling:

```
python paradys.py --patients P2 P3 --networks examples/networks_mini.csv --mutations examples/mutations_mini.csv --outdir results/ --kappa 2
```



## Output

In the output directory, we create the file `drivers.csv` storing detected significant driver genes for each analyzed patient. The column `'patient'` stores the ID of the patient that possesses the signifcant driver gene in `'driver'`, in combination with the associated dysregulation edge in `'dysregulation'` and the respective p-value of the chi-squared test in the column `'pvalue'`.

If additionally the `--scores` flag was set, the calculated impact scores for each driver are stored in the file `scores.csv`. The column `'patient'` stores the respective patient ID, `'driver'` the considered driver gene, and `'score'` the computed personalized impact score.



