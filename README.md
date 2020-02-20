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
N E X T F L O W  ~  version 20.01.0-edge
Launching `cosmo.nf` [ecstatic_leakey] - revision: b1e348ec8f
=========================================
COSMO => COrrection of Sample Mislabeling by Omics
=========================================
Usage:
nextflow run cosmo.nf
 Arguments:
   --pro_file              Protein expression data at gene level.
   --rna_file              RNA expressio data at gene level.
   --annotation_file       Sample annotation data.
   --annotation_attribute  Sample label(s) for prediction. Multiple labels
                           must be separated by ",".
   --method_id             1:SoonJye, 2:Sentieon, 3: 1+2. Default is 1.
   --task_id               The task ID, 2b or 2c, default is 2b.
   --out_dir               Output folder, default is "./output".
   --cpu                   The number of CPUs.
   --help                  Print help message.


```

### Input
The formats for both protein expression and mRNA expression data are the same. An example input of protein expression data (`--pro_file`) is shown below. The first column is the `gene ID` and all the other columns are the expression of proteins at gene level in different samples.

|  |Testing_1 | Testing_2 | Testing_3 | Testing_4 | Testing_5 | Testing_6 | Testing_7 | Testing_8 | Testing_9 | Testing_10 |        |
|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|------------|--------|--------|
| A1BG      | 1.5963    | 2.8484    | 2.1092    | 2.7922    | 2.4444    | 3.9907    | 3.6792    | 3.7321    | 3.6123     | 3.1739 |
| A2M       | 5.9429    | 5.0089    | 6.0823    | 6.0093    | 6.4553    | 6.0097    | 6.014     | 6.9721    | 4.4766     | 6.481  |
| AAAS      | 1.9337    | 2.951     | 3.5984    | 2.0419    | 2.1217    | 0.9662    | 1.0086    | NA        | 2.4936     | 2.2399 |
| AACS      | 1.7549    | NA        | 2.3948    | NA        | 0.9946    | 2.5969    | NA        | NA        | 1.6488     | NA     |
| AAGAB     | NA        | NA        | 0.9982    | NA        | 1.0282    | 1.6296    | NA        | NA        | 1.8141     | NA     |
| AAK1      | 1.0459    | 2.5435    | 1.7449    | NA        | 1.0653    | 0.9855    | 2.0395    | 1.1588    | NA         | NA     |


An example input of mRNA expression data (`--rna_file`) is shown below. The first column is the `gene ID` and all the other columns are the expression of proteins at gene level in different samples.

| |Testing_1 | Testing_2 | Testing_3 | Testing_4 | Testing_5 | Testing_6 | Testing_7 | Testing_8 | Testing_9 | Testing_10 |        |
|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|------------|--------|--------|
| A1BG      | 0.9394    | 1.0705    | 0.8846    | 1.5374    | 0.5717    | 1.3097    | 1.18      | 1.5178    | 0.1149     | 0.6317 |
| A1BG-AS1  | 0.2308    | 0.3563    | 0.1353    | 0.2406    | NA        | 0.4888    | 0.8254    | 0.1141    | 0.6046     | 0.0892 |
| A1CF      | 3.0512    | 2.9937    | 3.899     | 2.6612    | 3.0792    | 3.0881    | 2.1438    | 2.5478    | 2.9139     | 3.512  |
| A2M       | 6.2435    | 7.0947    | 6.232     | 5.3562    | 5.1192    | 7.3476    | 7.7089    | 7.3525    | 5.8959     | 5.8222 |
| A2M-AS1   | 0.9291    | 1.5417    | 0.399     | 0.8889    | 0.6882    | 1.7131    | 1.6557    | 1.9314    | 1.0677     | 1.3866 |
| A4GALT    | 1.9879    | 2.4708    | 1.4499    | 1.6679    | 1.8073    | 2.723     | 3.1722    | 3.5627    | 1.6291     | 2.1256 |

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
nextflow run cosmo.nf --pro_file example_data/for_2b/Testing_1/test_pro.tsv \
    --rna_file example_data/for_2b/Testing_1/test_rna.tsv \
    --sample_file example_data/for_2b/Testing_1/test_cli.tsv \
    --method_id 1 \
    --task_id 2b \
    --out_dir output
```
The data to run the above example can be found in this folder: "``example_data``".

### Output

