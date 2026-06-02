#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=5:00:00
#SBATCH --mem=50GB
#SBATCH --job-name=nextflow_rnaseq
#SBATCH --output=slurm_%j.out
#SBATCH --account=torch_pr_44_general

module purge

module load nextflow/26.04.0
	    
input=""
work=""
out=""

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -input) input="$2"; shift ;;
        -work) work="$2"; shift ;;
        -out) out="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

nextflow run main.nf --input "${input}" \
-w "${work}" \
--out_dir "${out}" \
--aligner "salmon"
	    


module purge