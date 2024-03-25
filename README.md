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

Follow the download instructions in the directory `examples/` to download exemplary preprocessed input files. You can then apply PARADYS on this input by calling:

```
python paradys.py --patients TCGA-CH-5766 --networks examples/networks.csv --mutations examples/mutations.csv --outdir results/ --kappa 3 --d 0.85 
```



## Output

In the output directory, we create a driver-dysregulation file for each input patient named `'<PATIENT_ID>_scores.csv'`. It consists of detected significant driver genes in the column `'driver'`, in combination with the associated dysregulation edge in `'dysregulation'` and the respective p-value of the chi-squared test in `'p-value'` .

If additionally the `--scores` flag was set, the calculated impact scores for each driver are stored in the patient-specific file `<PATIENT_ID>_scores.csv`.



