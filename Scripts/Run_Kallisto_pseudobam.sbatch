#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=5:00:00
#SBATCH --mem=32GB
#SBATCH --job-name=Kallsto_and_pseudobam
#SBATCH --output=slurm_%j.out
#SBATCH --account=pr_284_general

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
BASE=/path/to/base/analysis/directory/
OUT=${BASE}/rnaseq_analysis/${sample}
mkdir -p "$OUT"
cd "$BASE"

#purge all modules and load kallisto
module purge
module load fastqc/0.11.9
module load kallisto/0.46.1
module load trimgalore/0.6.6

#find if the analysis was single or paired-end. Names of the fastq files may need to be change case-by-case.
S1="${sample}.fastq.gz"
S1_paired="${sample}_1.fastq.gz"
S2="${sample}_2.fastq.gz"

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
    trim_galore --paired --fastqc --output_dir "$OUT" "$S1_file" "$S2_file"
    TRIM_S1="${OUT}/${sample}_1_val_1.fq.gz"
    TRIM_S2="${OUT}/${sample}_2_val_2.fq.gz"
else
    trim_galore --fastqc --output_dir "$OUT" "$S1_file"
    TRIM_S1="${OUT}/${sample}_trimmed.fq.gz"
fi

#run the kallisto quant; for single-end, need to know fragment length and sd of fragment length; default of both are 180 and 20
if [[ "$MODE" == "PE" ]]; then
    kallisto quant -i /scratch/projects/truonglab/GRCh38/${index} \
    -o "$OUT" \
    --pseudobam \
    "$TRIM_S1" "$TRIM_S2"
else
    kallisto quant -i /scratch/projects/truonglab/GRCh38/${index} \
    -o "$OUT" \
    --pseudobam \
    --single -l 180 -s 20 \
    "$TRIM_S1"
fi

module purge
