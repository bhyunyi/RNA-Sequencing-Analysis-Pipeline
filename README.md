# RNASeq Analysis Nextflow Pipeline

A reproducible Nextflow pipeline for bulk-RNASeq alignment and transcript count quantification. This workflow includes QC, trimming, alignment/pseudoalignment, quantification, and downstream count aggregation. This pipeline was developed in collaboration with David Truong's Lab at NYU. 
---

## Overview

This pipeline is for bioinformaticians to perform RNASeq alignment of multiple samples with a single command-line execution of a Nextflow workflow. Each sample is allowed a specification of a unique reference or index based on the transgene insertions or gene knockouts that are expected to be present in the sample. 

**Features**

- FASTQ quality control using FASTQC
- Adapter trimming using TRIMGALORE
- Transcript quantification using Salmon, STAR, Kallisto, or a combination of STAR and Salmon
- Multi-sample processing with the option to use unique indices or references for each sample
- SLURM/HPC support
- Containerized reproducibility (Singularity)
- Manual scripts for alignment that do not utilize Nextflow

---

## Workflow Overview

FastQC  
  ↓  
Trim Galore  
  ↓  
Salmon / STAR / Kallisto / STAR-Salmon  
  ↓  
Quantification  

---

## Repository Structure

```text
repo/
├── Nextflow/
│   ├── main.nf
│   ├── nextflow.config
│   ├── run_RNASeq_nextflow.sbatch
│   └── samplesheet.csv
│
├── manual_scripts/
│   ├── Kallisto_Align
│   └── STAR_Align
│
└── README.md
```

## Usage Guide

### 1. Download and Edit Scripts

- main.nf: Includes each process and the workflow to execute the processes. Does not need to be changed from the original version in this repo.
- nextflow.config: includes specifications for input parameters, singularity parameters, process defaults, and execution specifications. Some lines must be changed from the original version in order to be used.
  - All paths to reference files and indices must be changed to the user's respective paths.
  - Paths to the user's saved singularity overlay and container image must be changed to the user's respecitve paths. 
- samplesheet.csv: includes a comma separated list in which each row represents an RNASeq sample to be aligned and analyzed. Each row has four columns that must be filled in based on the samples to be analyzed.
  - Sample Sheet Column Descriptions:

| Column | Description |
|---|---|
| `sample1` | Path to first FASTQ file. Only fill in this column if experiment is single-end |
| `sample2` | Path to second FASTQ file. Leave this column blank if experiment is single-end |
| `index` | Name of the index specified in the params of the config file. |
| `vendor` | Vendor in which the RNASeq was performed. For now, there are only two options of 'novogene' or 'plasmidsaurus' |

- run_RNASeq_nextflow.sbatch: includes the sbatch script to run nextflow given that the main.nf, nextflow.config, and samplesheet.csv are in the same directory. This is where the aligner specification is made. 

### 2. Run Sbatch Script

Run the run_RNASeq_nextflow.sbatch on the command line with paths for the input samplesheet, nextflow work directory, and output directory. Ensure that the main.nf and nextflow.config files are within the same directory in which this sbatch file is being called. The work and output directories do not need to exist for the command to execute. 

```bash
sbatch run_RNASeq_nextflow.sbatch -input /path/to/samplesheet -work /path/to/nextflow_work/ -out /path/to/nextflow_out/
```

### 3. Check Outputs

Within the path that was specified for the nextflow output directory, you will find the following structure:
```text
nextflow_output/
├── sample/
│   ├── fastqc/
│   ├── aligner/
│   └── trimmed
```
The fastqc folder will contain the QC html reports for the FASTQ file(s) for the sample. The aligner folder will be named kallisto, salmon, star, or star-salmon and will contain the transcript count file and other outputs from the specific aligner. The trimmed folder will contain the trimmed FASTQ file(s). 

### 4. Downstream Differential Expression Analysis with Streamlit App

The transcript quantification outputs from this pipeline are meant to be used with the streamlit app shown in the following repo: https://github.com/am15443/rnaseq-app

The user may have the option to perform their own Differential Expression Analysis using the R Scripts for DESeq2 in the manual_scripts folder of this repo. Currently, there are only manual scripts from Kallisto and STAR aligners.
