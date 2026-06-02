#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def helpMessage() {
    if (params.help) {
        log.info """
        =========================================
        RNA-seq Alignment Pipeline
        =========================================
        Usage:
            nextflow run main.nf \
            --input samplesheet.csv \
            --out_dir /path/to/output/directory/ \
            --aligner "star"
    
        Samplesheet CSV columns:
            sample1,sample2,index,vendor
    
        Options:
            --input      Path to samplesheet CSV  [required]
            --out_dir    Output directory
            --aligner    Aligner to be used
            --help       Show this help message
            
        Aligner Options:
           "star"        Default aligner
           "salmon"   
           "kallisto"
           "star_salmon"
        =========================================
        """.stripIndent()
    }
}


// Processes ==================================================


process FASTQC {
    tag "${sample}"
    publishDir params.out_dir, mode: 'copy', saveAs: { filename -> "${sample}/fastqc/${filename}" }

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("*.html")                                                                      , emit: html
    tuple val(sample), path("*.zip")                                                                       , emit: zip
    path "fastqc_version.txt"                                                                              , emit: fastqc_version

    script:
    def read_args = reads instanceof List ? reads.join(' ') : reads
    """
    source /ext3/env.sh
    fastqc ${read_args} -o .
    
    fastqc --version > fastqc_version.txt

    """
}

process TRIM_GALORE {
    tag "${sample}"
    publishDir params.out_dir, mode: 'copy', saveAs: { filename -> "${sample}/trimmed/${filename}" }

    input:
    tuple val(sample), val(mode), path(reads)

    output:
    tuple val(sample), val(mode), path("*.{fq,fq.gz,fastq,fastq.gz}")                                              , emit: trimmed_reads
    tuple val(sample), path("*_trimming_report.txt")                                                               , emit: trim_reports
    path "trim_galore_version.txt"                                                                                 , emit: trim_galore_version


    script:
    if (mode == "PE") {
        def r1 = reads[0]
        def r2 = reads[1]
        """
        source /ext3/env.sh
        conda run -n trim_env trim_galore --paired --fastqc --output_dir . ${r1} ${r2}
        
        conda run -n trim_env trim_galore --version > trim_galore_version.txt
        """
    } else {
        """
        source /ext3/env.sh
        conda run -n trim_env trim_galore --fastqc --output_dir . ${reads}
        
        conda run -n trim_env trim_galore --version > trim_galore_version.txt

        """
    }
}

process STAR_ALIGN {
    tag "${sample}"
    publishDir params.out_dir, mode: 'copy', saveAs: { filename -> "${sample}/star/${filename}" }

    input:
    tuple val(sample), val(mode), path(trimmed_reads), val(index)

    output:
    tuple val(sample), val(index), path("Aligned.sortedByCoord.out.bam")     , emit: bam
    tuple val(sample), path("ReadsPerGene.out.tab")                          , emit: counts
    tuple val(sample), path("Log.final.out")                                 , emit: log
    path "STAR_version.txt"                                                  , emit: STAR_version

    script:
    def genome_dir = params.references[index]['star_index']
    def reads_arg  = (mode == "PE") ? "${trimmed_reads[0]} ${trimmed_reads[1]}" : "${trimmed_reads}"
    def zcat_arg   = reads_arg.contains('.gz') ? "--readFilesCommand zcat" : ""
    """
    source /ext3/env.sh
    STAR --runThreadN ${task.cpus} \\
         --readFilesIn ${reads_arg} \\
         ${zcat_arg} \\
         --genomeDir ${genome_dir} \\
         --quantMode GeneCounts \\
         --outSAMtype BAM SortedByCoordinate \\
         --outFileNamePrefix ./
         
    STAR --version > STAR_version.txt
    """
}

process STAR_ALIGN_TRANSCRIPTOME {
    tag "${sample}"
    publishDir params.out_dir, mode: 'copy', saveAs: { filename -> "${sample}/star/${filename}" }
    
    input:
    tuple val(sample), val(mode), path(trimmed_reads), val(index)
    
    output:
    tuple val(sample), val(index), path("Aligned.toTranscriptome.out.bam")   , emit: bam
    tuple val(sample), path("Log.final.out")                                 , emit: log
    path "STAR_version.txt"                                                  , emit: STAR_version


    script:
    def genome_dir = params.references[index]['star_index']
    def reads_arg  = (mode == "PE") ? "${trimmed_reads[0]} ${trimmed_reads[1]}" : "${trimmed_reads}"
    def zcat_arg   = reads_arg.contains('.gz') ? "--readFilesCommand zcat" : ""
    """
    source /ext3/env.sh
    STAR --runThreadN ${task.cpus} \\
         --readFilesIn ${reads_arg} \\
         ${zcat_arg} \\
         --genomeDir ${genome_dir} \\
         --quantMode TranscriptomeSAM \\
         --outSAMtype None \\
         --outFileNamePrefix ./
         
    STAR --version > STAR_version.txt
    """
}

process SALMON_STAR_QUANT {
    tag "${sample}"
    publishDir params.out_dir, mode: 'copy', saveAs: { filename -> "${sample}/salmon/${filename}" }
    
    input:
    tuple val(sample), val(index), path(transcriptome_bam)
    
    output:
    tuple val(sample), path("${sample}_salmon/quant.sf")                       , emit: quant
    path "Salmon_version.txt"                                                  , emit: Salmon_version

    script:
    def transc_fa = params.references[index]['transc_fa']

    """
    source /ext3/env.sh
    salmon quant \
            -t ${transc_fa} \
            -l A \
            -a ${transcriptome_bam} \
            -p ${task.cpus} \
            -o ${sample}_salmon
        
    salmon --version > Salmon_version.txt
    """
    
}

process SALMON_ALIGN {
    tag "${sample}"
    publishDir params.out_dir, mode: 'copy', saveAs: { filename -> "${sample}/salmon/${filename}" }
    
    input:
    tuple val(sample), val(mode), path(trimmed_reads), val(index)
    
    output:
    tuple val(sample), val(index), path("quant.sf")                             , emit: sf
    tuple val(sample), val(index), path("logs")                                 , emit: logs
    tuple val(sample), val(index), path("aux_info")                             , emit: aux_info
    tuple val(sample), val(index), path("cmd_info.json")                        , emit: cmd_json
    tuple val(sample), val(index), path("lib_format_counts.json")               , emit: counts_json
    path "Salmon_version.txt"                                                   , emit: Salmon_version


    script:
    def index_path = params.references[index]['salmon_index']
    def reads_arg  = (mode == "PE") ? "-1 ${trimmed_reads[0]} -2 ${trimmed_reads[1]}" : "-r ${trimmed_reads}"
    """
    source /ext3/env.sh
    salmon quant -i ${index_path} \\
        -l A \\
        ${reads_arg} \\
        -p ${task.cpus} \\
        --validateMappings \\
        --gcBias \\
        --seqBias \\
        -o .
        
    salmon --version > Salmon_version.txt
    """
}

process KALLISTO_ALIGN{
    tag "${sample}"
    publishDir params.out_dir, mode: 'copy', saveAs: { filename -> "${sample}/kallisto/${filename}" }
    
    input:
    tuple val(sample), val(mode), path(trimmed_reads), val(index)
    
    output:
    tuple val(sample), val(index), path("abundance.h5")                                           , emit: h5
    tuple val(sample), val(index), path("abundance.tsv")                                          , emit: tsv
    tuple val(sample), val(index), path("pseudoalignments.bam")                                   , emit: bam
    tuple val(sample), val(index), path("run_info.json")                                          , emit: json
    path "Kallisto_version.txt"                                                                   , emit: Kallisto_version

    script:
    def index_path = params.references[index]['kallisto_index']
    def reads_arg  = (mode == "PE") ? "${trimmed_reads[0]} ${trimmed_reads[1]}" : "${trimmed_reads}"
    def frag_args  = (mode == "SE") ? "--single -l 250 -s 25" : ""
    """
    source /ext3/env.sh
    conda run -n kallisto_0461 kallisto quant \\
        -i ${index_path} \\
        -o . \\
        --pseudobam \\
        ${frag_args} \\
        ${reads_arg}
        
    conda run -n kallisto_0461 kallisto version > Kallisto_version.txt     
    """
}


process UMI_DEDUP {
    tag "${sample}"
    publishDir params.out_dir, mode: 'copy', saveAs: { filename -> "${sample}/dedup/${filename}" }

    input:
    tuple val(sample), val(index), path(bam), val(vendor)

    output:
    tuple val(sample), val(index), path("dedup.bam"), val(vendor)                                                 , emit: bam
    path "UMI_collapse_version.txt"                                                                               , emit: UMI_collapse_version


    script:
    """
    source /ext3/env.sh
    java -Xmx48G -jar /ext3/miniforge3/share/umicollapse-1.1.0-0/umicollapse.jar bam \\
        -i ${bam} \\
        -o dedup.bam \\
        --umi-sep _
        
    echo "umicollapse: 1.1.0" > UMI_collapse_version.txt
    """
}

process FEATURE_COUNTS {
    tag "${sample}"
    publishDir params.out_dir, mode: 'copy', saveAs: { filename -> "${sample}/counts/${filename}" }
    cpus 4

    input:
    tuple val(sample), val(index), path(bam), val(vendor), val(mode)

    output:
    tuple val(sample), path("gene_feature_counts.txt")                                                                   , emit: counts
    tuple val(sample), path("gene_feature_counts.txt.summary")                                                           , emit: summary
    path "Feature_counts_version.txt"                                                                                    , emit: Feature_counts_version


    script:
    def gtf
    if (index.contains("PGP1")) {
        gtf = params.references['PGP1_genome']['gtf']
    } else if (index.contains("miniNfull")) {
        gtf = params.references['miniNfull_genome']['gtf']
    } else {
        error "Unknown index type '${index}' for sample ${sample}. Cannot determine GTF."
    }
    def end_args = (mode == "PE") ? "-p -B -C" : ""

    """
    source /ext3/env.sh
    featureCounts \\
        -T ${task.cpus} \\
        -t exon \\
        -g gene_id \\
        -a ${gtf} \\
        -o gene_feature_counts.txt \\
        ${end_args} \\
        ${bam}
        
    featureCounts -v 2> Feature_counts_version.txt
    """
}


// Helper function to parse FASTq paths

def resolve_reads(sample1, sample2) { 
    if (sample2 == null || sample2.toString().trim() == "") {
        def se = file(sample1)
        if (se.exists()) return ["SE", se]
        else error "No FASTQ file found for SE sample: '${sample1}'"
    }

    def r1 = file(sample1)
    def r2 = file(sample2)
    if (r1.exists() && r2.exists()) return ["PE", [r1, r2]]

    error "No FASTQ files found for PE sample: '${sample1}' / '${sample2}'"
}


// Workflow ==================================================

workflow {

    if (params.help) {
        helpMessage()
        exit 0
    }

    if (!params.input)    error "ERROR: --input is required"
    if (!params.out_dir)  error "ERROR: --out_dir is required"

    // Parse samplesheet: sample,index,vendor
    Channel
        .fromPath(params.input)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def ref = params.references[row.index]
            if (!ref) error "Unknown reference: ${row.index}"
            def (mode, reads) = resolve_reads(row.sample1, row.sample2)
            // derive a clean ID from the filename
            def sample_id = file(row.sample1).simpleName.replaceAll(/(_1|_R1)$/, "")
            tuple(sample_id, mode, reads, row.index, row.vendor)
        }
        .set { samples_ch }
        
    reads_for_fastqc = samples_ch.map { sample_id, mode, reads, index, vendor ->
        tuple(sample_id, reads)
    }
    
    reads_for_trim = samples_ch.map { sample_id, mode, reads, index, vendor ->
        tuple(sample_id, mode, reads)
    }
    
    FASTQC(reads_for_fastqc)
    TRIM_GALORE(reads_for_trim)

    // Attach index back to trimmed reads for STAR
    trimmed_with_meta = TRIM_GALORE.out.trimmed_reads
        .join(
            samples_ch.map { sample_id, mode, reads, index, vendor -> tuple(sample_id, index) }
        )

    if (params.aligner == "star") {
        STAR_ALIGN(trimmed_with_meta)
    
        // Attach vendor to BAM for dedup decision
        bam_with_vendor = STAR_ALIGN.out.bam
            .join(
                samples_ch.map { sample_id, mode, reads, index, vendor -> tuple(sample_id, vendor) }
            )

        // split into two channels based on vendor
        plasmidsaurus_ch = bam_with_vendor.filter { sample_id, index, bam, vendor -> vendor == "plasmidsaurus" }
        novogene_ch      = bam_with_vendor.filter { sample_id, index, bam, vendor -> vendor == "novogene" }
        
        // only plasmidsaurus goes through dedup
        UMI_DEDUP(plasmidsaurus_ch)
        
        // merge both back together for featureCounts
        all_bams_ch = UMI_DEDUP.out.bam.mix(novogene_ch)   
        
        all_bams_with_mode = all_bams_ch
            .join(
                samples_ch.map { sample_id, mode, reads, index, vendor -> tuple(sample_id, mode) }
            )

        FEATURE_COUNTS(all_bams_with_mode)
        
    }
    
    if (params.aligner == "star_salmon") {
        STAR_ALIGN_TRANSCRIPTOME(trimmed_with_meta)
        
        transcriptome_bam = STAR_ALIGN_TRANSCRIPTOME.out.bam
        
        SALMON_STAR_QUANT(transcriptome_bam)
    }
    
    if (params.aligner == "salmon") {
        SALMON_ALIGN(trimmed_with_meta)
    }
    
    if (params.aligner == "kallisto") {
        KALLISTO_ALIGN(trimmed_with_meta)
    }
}