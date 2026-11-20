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

 - Sample Sheet example:

```text
sample1,sample2,index,vendor
/path/to/sample/fastq1.fq,/path/to/sample/fastq2.fq,index_name,vendor_name
```

- The 'index' field in the Sample Sheet refers to the name of the unique transgene index
that is listed within the "references" dictionary within "nextflow.config" file. The 
required indices and gtf  files for STAR, Kallisto and Salmon must be generated beforehand,
 and the paths to the relevant files must be listed with the respective keys in the 
"references" dictionary.
 - Below is and example of the "references" dictionary within the "nextflow.config" file
(PGP1_genome, miniNfull_genome, and GRCh38 are the identifiers to use for the 'index' field):

```text
    references = [
        'PGP1_genome': [
            star_index:     '/path/to/PGP1 Star index',
            kallisto_index: '/path/to/PGP1 kallisto index',
            salmon_index:   '/path/to/PGP1 salmon index',
            gtf:            '/path/to/PGP1 gtf file',
            transc_fa:      '/path/to/PGP1 transcriptome fasta'
        ],
        'miniNfull_genome': [
            star_index:     '/path/to/miniNfull Star index',
            kallisto_index: '/path/to/miniNfull kallisto index',
            salmon_index:   '/path/to/miniNfull salmon index',
            gtf:            '/path/to/miniNfull gtf file',
            transc_fa:      '/path/to/miniNfull transcriptome fasta'
        ],
        'GRCh38': [
            star_index:     '/path/to/GRCh38 Star index',
            kallisto_index: '/path/to/GRCh38 kallisto index',
            salmon_index:   '/path/to/GRCh38 salmon index',
            gtf:            '/path/to/GRCh38 gtf file',
            transc_fa:      '/path/to/GRCh38 transcriptome fasta'
        ]
        
    ]
``` 

- The 'vendor' field in the samplesheet can be either "plasmidsaurus" or "novogene" since 
these are the two vendors that the Truong Lab utilizes for bulk RNA-sequencing. 
 - plasmidsaurus: single-end read analysis with UMI deduplication (plasmidsaurus utilizes 
UMIs in their single-end, long read sequencing) 
 - novogene: paired-end read analysis without UMI depublication

- run_RNASeq_nextflow.sbatch: includes the sbatch script to run nextflow given that the 
main.nf, nextflow.config, and samplesheet.csv are in the same directory. This is where 
the aligner specification is made.

### 2. Run Sbatch Script

Run the run_RNASeq_nextflow.sbatch on the command line with paths for the input 
samplesheet, nextflow work directory, and output directory. Ensure that the main.nf 
and nextflow.config files are within the same directory in which this sbatch file is 
being called. The work and output directories do not need to exist for the command to 
execute. 

```bash
sbatch run_RNASeq_nextflow.sbatch -input /path/to/samplesheet -work /path/to/nextflow_work/ -out /path/to/nextflow_out/
```

To switch aligners, you must go directly into the run_RNASeq_nextflow.sbatch script and 
edit the designation for the "--aligner" command to "salmon", "star", or "kallisto".

```text
nextflow run main.nf --input "${input}" \
-w "${work}" \
--out_dir "${out}" \
--aligner "salmon or star or kallisto"
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

### 4. HPC and Slurm Allocations For Each Process Within the Pipeline

- Below are the SLURM Scheduler allocations for the Nextflow Pipeline:

```bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=5:00:00
#SBATCH --mem=50GB
#SBATCH --job-name=nextflow_rnaseq
#SBATCH --output=slurm_%j.out
```

- The number of CPUs and RAM allocated for the process is in line with the CPU and RAM 
allocations for each part of the pipeline that is within the "nextflow.config".

```text
// ── Process Defaults
process {
    executor       = 'slurm'
    clusterOptions = '--account=torch_pr_44_general'
    container      = '/share/apps/images/cuda12.1.1-cudnn8.9.0-devel-ubuntu22.04.2.sif'

    // Defaults for all processes
    cpus   = 1
    memory = '8 GB'
    time   = '1h'

    // Per-process resource overrides
    withName: 'FASTQC' {
        cpus   = 2
        memory = '8 GB'
        time   = '1h'
    }
    withName: 'TRIM_GALORE' {
        cpus   = 2
        memory = '30 GB'
        time   = '3h'
    }
    withName: 'STAR_ALIGN' {
        cpus   = 4
        memory = '40 GB'
        time   = '3h'
    }
    withName: 'UMI_DEDUP' {
        cpus   = 4
        memory = '45 GB' 
        time   = '2h'
    }
    withName: 'FEATURE_COUNTS' {
        cpus   = 4
        memory = '10 GB'
        time   = '1h'
    }
    withName: 'STAR_ALIGN_TRANSCRIPTOME' {
        cpus   = 4
        memory = '40 GB'
        time   = '3h'
    }
    withName: 'KALLISTO_ALIGN' {
        cpus   = 4
        memory = '12 GB'
        time   = '2h'
    }
    withName: 'SALMON_STAR_QUANT' {
        cpus   = 4
        memory = '12 GB'
        time   = '1h'
    }
    withName: 'SALMON_ALIGN' {
    cpus   = 4
    memory = '15 GB'
    time   = '2h'
    }   
}

// ── Execution Settings 
executor {
    name           = 'slurm'
    queueSize      = 50     // max concurrent jobs
    submitRateLimit = '10/1min'
}

``` 

- No task has more than 50GB of RAM or 4 CPUs allocated. Nextflow is limited to start 
50 jobs at a time on the SLURM scheduler with a limit of submitting 10 jobs per minute. 

### 5. Downstream Differential Expression Analysis with Streamlit App

The transcript quantification outputs from this pipeline are meant to be used with the streamlit app shown in the following repo: https://github.com/am15443/rnaseq-app

The user may have the option to perform their own Differential Expression Analysis using the R Scripts for DESeq2 in the manual_scripts folder of this repo. Currently, there are only manual scripts from Kallisto and STAR aligners.

### 6. Downstream Differential Expression Analysis with /manual_scripts/


