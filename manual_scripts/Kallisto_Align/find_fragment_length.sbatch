#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --time=5:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=Salmon_frag_length
#SBATCH --output=slurm_%j.out
#SBATCH --account=torch_pr_44_general


module purge

singularity exec --nv \
	    --overlay /path/to/environment/envir.ext3:ro \
	    /share/apps/images/cuda12.1.1-cudnn8.9.0-devel-ubuntu22.04.2.sif \
	    /bin/bash -c 'source /ext3/env.sh

#initialize variables
sample=""

#parse command-line arguments for -idx and -sra flags
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -sample) sample="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

cd /path/to/frag_length/output/

#run the salmon pseudoalignment
salmon quant -i /projects/rps/truond01/truonglab/GRCh38/Salmon_idx \
-l A \
-r /path/to/Raw_data/${sample}.fastq.gz \
--validateMappings -o salmon_tmp

python3 <<EOF
import gzip
import numpy as np
import os

fld_path = "salmon_tmp/aux_info/fld.gz"

if not os.path.exists(fld_path):
    print("ERROR: fld.gz not found — fragment inference may have failed.")
    exit(1)

with gzip.open(fld_path, "rb") as f:
    counts = np.frombuffer(f.read(), dtype=np.int32)

total = counts.sum()

if total == 0:
    print("ERROR: No fragment counts detected.")
    exit(1)

lengths = np.arange(len(counts))

mean = (lengths * counts).sum() / total
sd = np.sqrt((counts * (lengths - mean)**2).sum() / total)

print(f"Fragment length mean: {mean:.2f}")
print(f"Fragment length SD:   {sd:.2f}")
EOF

' bash "$@"

module purge
