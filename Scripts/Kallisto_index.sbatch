#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=5:00:00
#SBATCH --mem=32GB
#SBATCH --job-name=Kallsto_index
#SBATCH --output=slurm_%j.out
#SBATCH --account=torch_pr_44_general

module purge

singularity exec --nv \
	    --overlay /scratch/by2372/ML_envir/ML_envir.ext3:ro \
	    /share/apps/images/cuda12.1.1-cudnn8.9.0-devel-ubuntu22.04.2.sif \
	    /bin/bash -c 'source /ext3/env.sh

cd /scratch/by2372/Microglia_RNA_Seq_10FEB2026/References/

conda run -n kallisto_0461 kallisto index -i New.GRCh38.cdna.all.idx New.GRCh38.cdna.all.fa

' bash "$@"

module purge