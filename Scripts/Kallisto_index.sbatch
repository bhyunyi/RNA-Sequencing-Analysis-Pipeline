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
	    --overlay /path/to/environment/ML_envir.ext3:ro \
	    /share/apps/images/cuda12.1.1-cudnn8.9.0-devel-ubuntu22.04.2.sif \
	    /bin/bash -c 'source /ext3/env.sh

cd /path/to/References/

conda run -n kallisto_0461 kallisto index -i New.idx reference.fa

' bash "$@"

module purge
