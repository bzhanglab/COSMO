# COSMO

## Overview

**COSMO: COrrection of Sample Mislabeling by Omics**

Multi-omics Enabled Sample Mislabeling Correction

## Installation

1. Download neoflow:

```sh
git clone https://github.com/bzhanglab/cosmo
```

2. Install [Docker](https://docs.docker.com/install/) (>=19.03).

3. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html). More information can be found in the Nextflow [get started](https://www.nextflow.io/docs/latest/getstarted.html) page.

All other tools used by COSMO have been dockerized and will be automatically installed when COSMO is run in the first time on a computer.

## Usage

```sh
○ → nextflow run cosmo.nf --help
N E X T F L O W  ~  version 21.02.0-edge
Launching `cosmo.nf` [high_noyce] - revision: 46dc5c6c96
=========================================
COSMO => COrrection of Sample Mislabeling by Omics
=========================================
Usage:
nextflow run cosmo.nf
 Arguments:
   --d1_file               Dataset with quantification data at gene level.
   --d2_file               Dataset with quantification data at gene level.
   --d1_type               The type for dataset d1. This is used to label the dataset in the output files.
   --d2_type               The type for dataset d2. This is used to label the dataset in the output files.
   --cli_file              Sample annotation data.
   --cli_attribute         Sample attribute(s) for prediction. Multiple attributes
                           must be separated by ",".
   --out_dir               Output folder, default is "./output".
   --threads               The number of threads.
   --help                  Print help message.


```

### Input
The formats for both datasets (`--d1_file`, `--d2_file`) are the same. An example input of quantification dataset (`--d1_file` or `--d2_file`) is shown below. The first column is the `gene ID` and all the other columns are the expression of proteins at gene level in different samples.

|  |Testing_1 | Testing_2 | Testing_3 | Testing_4 | Testing_5 | Testing_6 | Testing_7 | Testing_8 | Testing_9 | Testing_10 |        |
|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|------------|--------|--------|
| A1BG      | 1.5963    | 2.8484    | 2.1092    | 2.7922    | 2.4444    | 3.9907    | 3.6792    | 3.7321    | 3.6123     | 3.1739 |
| A2M       | 5.9429    | 5.0089    | 6.0823    | 6.0093    | 6.4553    | 6.0097    | 6.014     | 6.9721    | 4.4766     | 6.481  |
| AAAS      | 1.9337    | 2.951     | 3.5984    | 2.0419    | 2.1217    | 0.9662    | 1.0086    | NA        | 2.4936     | 2.2399 |
| AACS      | 1.7549    | NA        | 2.3948    | NA        | 0.9946    | 2.5969    | NA        | NA        | 1.6488     | NA     |
| AAGAB     | NA        | NA        | 0.9982    | NA        | 1.0282    | 1.6296    | NA        | NA        | 1.8141     | NA     |
| AAK1      | 1.0459    | 2.5435    | 1.7449    | NA        | 1.0653    | 0.9855    | 2.0395    | 1.1588    | NA         | NA     |


The input for parameter `--sample_file` is the sample annotation file and an example is shown below:


| sample    | age | gender | stage | colon_rectum | msi         | tumor_normal |
|-----------|-----|--------|-------|--------------|-------------|--------------|
| Testing_1 | 47  | Female | High  | Colon        | MSI-Low/MSS | Tumor        |
| Testing_2 | 68  | Female | High  | Rectum       | MSI-Low/MSS | Tumor        |
| Testing_3 | 52  | Male   | Low   | Colon        | MSI-Low/MSS | Tumor        |
| Testing_4 | 54  | Female | Low   | Colon        | MSI-High    | Tumor        |
| Testing_5 | 72  | Male   | High  | Colon        | MSI-Low/MSS | Tumor        |
| Testing_6 | 61  | Male   | High  | Colon        | MSI-Low/MSS | Tumor        |
| Testing_7 | 58  | Female | High  | Colon        | MSI-High    | Tumor        |
| Testing_8 | 73  | Male   | Low   | Colon        | MSI-Low/MSS | Tumor        |
| Testing_9 | 68  | Male   | Low   | Colon        | MSI-Low/MSS | Tumor        |


Below is an example run COSMO:
```sh
nextflow run cosmo.nf --d1_file example_data/test_pro.tsv \
    --d2_file example_data/test_rna.tsv \
    --cli_file example_data/test_cli.tsv \
    --cli_attribute "gender,msi" \
    --out_dir output
```
The data to run the above example can be found in this folder: "``example_data``".

### Output

