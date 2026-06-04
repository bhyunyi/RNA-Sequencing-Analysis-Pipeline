# RNASeq Analysis Nextflow Pipeline

A reproducible Nextflow pipeline for bulk-RNASeq alignment and transcript count quantification. This workflow includes QC, trimming, alignment/pseudoalignment, quantification, and downstream count aggregation.
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

---

## Workflow Overview

FASTQ
  ↓
FastQC
  ↓
Trim Galore
  ↓
Salmon / STAR / Kallisto / STAR-Salmon
  ↓
Quantification
  ↓
Count Matrix Generation

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

