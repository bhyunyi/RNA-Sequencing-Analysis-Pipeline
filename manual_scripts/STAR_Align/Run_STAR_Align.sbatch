#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=5:00:00
#SBATCH --mem=55GB
#SBATCH --job-name=STAR_aln
#SBATCH --output=slurm_%j.out
#SBATCH --account=torch_pr_44_general

module purge

singularity exec --nv \
	    --overlay /scratch/by2372/ML_envir/ML_envir.ext3:ro \
	    /share/apps/images/cuda12.1.1-cudnn8.9.0-devel-ubuntu22.04.2.sif \
	    /bin/bash -c 'source /ext3/env.sh

#initialize variables

sample=""
index=""
vendor=""

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -sample) sample="$2"; shift ;;
        -index) index="$2"; shift ;;
        -vendor) vendor="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

echo "Running sample=$sample index=$index vendor=$vendor"

#define the paths to be used and make directories where needed
BASE=/projects/rps/truond01/truonglab/GCBW8H_Flp_out_p55_iFNg/Raw_data/
OUT=/scratch/by2372/Flp_out_lines_STAR_outputs/${sample}
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
else
    echo "No FASTQ files found for sample $sample" >&2
    exit 1
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

if [[ "$MODE" == "PE" ]]; then
    STAR --runThreadN 4 \
    --readFilesIn "$TRIM_S1" "$TRIM_S2" \
    --genomeDir /scratch/by2372/FlpOut_lines_p55_RNASEQ_10MAR2026/STAR_ref/${index}/ \
    --quantMode GeneCounts \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${OUT}/
else
    STAR --runThreadN 4 \
    --readFilesIn "$TRIM_S1" \
    --genomeDir /scratch/by2372/FlpOut_lines_p55_RNASEQ_10MAR2026/STAR_ref/${index}/ \
    --quantMode GeneCounts \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${OUT}/
fi

cd "$OUT"

BAM=""

if [[ "$vendor" == "plasmidsaurus" ]]; then
    java -Xmx48G -jar /ext3/miniforge3/share/umicollapse-1.1.0-0/umicollapse.jar bam \
      -i Aligned.sortedByCoord.out.bam \
      -o dedup.bam \
      --umi-sep _
      
    BAM="dedup.bam"
elif [[ "$vendor" == "novogene" ]]; then
    BAM="Aligned.sortedByCoord.out.bam"
else
    echo "Invalid Vendor Name!!" >&2
    exit 1
fi
  
if [[ "$index" == *PGP1* ]]; then
    GTF="GRCh38_reference_PGP1.gtf"
elif [[ "$index" == *miniNfull* ]]; then
    GTF="GRCh38_reference_miniNfull.gtf"
else
    echo "Unknown index type: $index" >&2
    exit 1
fi
  
featureCounts \
  -T 4 \
  -t exon \
  -g gene_id \
  -a /scratch/by2372/FlpOut_lines_p55_RNASEQ_10MAR2026/STAR_ref/${GTF} \
  -o gene_feature_counts.txt \
  "$BAM"

' bash "$@"

module purge
