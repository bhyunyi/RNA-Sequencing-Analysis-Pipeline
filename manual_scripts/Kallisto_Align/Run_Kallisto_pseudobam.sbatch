#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=5:00:00
#SBATCH --mem=32GB
#SBATCH --job-name=Kallsto_and_pseudobam
#SBATCH --output=slurm_%j.out
#SBATCH --account=pr_284_general

module purge

singularity exec --nv \
	    --overlay /path/to/environment/envir.ext3:ro \
	    /share/apps/images/cuda12.1.1-cudnn8.9.0-devel-ubuntu22.04.2.sif \
	    /bin/bash -c 'source /ext3/env.sh

#initialize variables

sample=""
index="Homo_sapiens.GRCh38.cdna.all_kallisto.idx"

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -sample) sample="$2"; shift ;;
        -index) index="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

#define the paths to be used and make directories where needed
BASE=/path/to/sample/fastqs/${sample}
OUT=/path/to/output/directory/${sample}
mkdir -p "$OUT"
cd "$BASE"

#find if the analysis was single or paired-end. Names of the fastq files may need to be change case-by-case.
S1="${sample}.fastq"
S1_paired="${sample}_1.fastq"
S2="${sample}_2.fastq"

if [[ -f "$S1_paired" && -f "$S2" ]]; then
    MODE="PE"
    S1_file="$S1_paired"
    S2_file="$S2"
elif [[ -f "$S1" ]]; then
    MODE="SE"
    S1_file="$S1"
fi

#fastqc the fastq.
if [[ "$MODE" == "PE" ]]; then
    fastqc "$S1_file" "$S2_file" -o "$OUT"
else
    fastqc "$S1_file" -o "$OUT"
fi

#trim the reads using trim galore
if [[ "$MODE" == "PE" ]]; then
    conda run -n trim_env trim_galore --paired --fastqc --output_dir "$OUT" "$S1_file" "$S2_file"
    TRIM_S1="${OUT}/${sample}_1_val_1.fq"
    TRIM_S2="${OUT}/${sample}_2_val_2.fq"
else
    conda run -n trim_env trim_galore --fastqc --output_dir "$OUT" "$S1_file"
    TRIM_S1="${OUT}/${sample}_trimmed.fq"
fi

#run the kallisto quant; for single-end, need to know fragment length and sd of fragment length; default of both are 180 and 20
if [[ "$MODE" == "PE" ]]; then
    conda run -n kallisto_0461 kallisto quant -i /scratch/projects/truonglab/GRCh38/${index} \
    -o "$OUT" \
    --pseudobam \
    "$TRIM_S1" "$TRIM_S2"
else
    conda run -n kallisto_0461 kallisto quant -i /scratch/projects/truonglab/GRCh38/${index} \
    -o "$OUT" \
    --pseudobam \
    --single -l 180 -s 20 \
    "$TRIM_S1"
fi

' bash "$@"

module purge
