# Below are the directions to use this RNA-Seq analysis Pipeline.

## Part 1: Running Kallisto Pseudoalignment and Pseudobam generation

The “Run_Kallisto_pseudobam.sbatch” takes an input of an RNA sequenced sample name, and an optional input of a kallisto reference index name. First, FastQC is used to perform quality control on the fastq files of the sample. Then, Trimgalore is used to trim the reads of the raw data. Trimgalore performs another round of FastQC on the trimmed reads and outputs the trimmed reads into a specified output folder. These trimmed reads are analyzed using Kallisto, resulting in outputs that specify transcript counts and a pseudobam file. 

Change the following lines in “Run_Kallisto_pseudobam.sbatch” based on the RNA sequencing experiment and location of the raw data and kallisto index. 

Line 14: index="Homo_sapiens.GRCh38.cdna.all_kallisto.idx"
- The default index is set to a kallisto index of the human genome GRCh38 that lives in the truong lab projects directory. Lines 70 and 75 call on this kallisto index in which the location in the truong lab projects folder is specified. Change this index file name to the desired index if the default human genome index is not meant to be used for a specific RNA sequencing experiment. Change the index directory locations in lines 70 and 75 if the locations are changed. You can also specify a different index name within the same directory using the “-index” call at the command line when executing this sbatch script. 

Line 26: BASE=/path/to/sample/fastqs/${sample}
- Input the path to the raw data fastq files where it says “/path/to/sample/fastqs/”. Keep the “${sample}” variable as is at the end of the file path. 

Line 27: OUT=/path/to/output/directory/${sample}
- Input the path in which the RNA seq outputs will be deposited after analysis is performed. Replace “/path/to/output/directory/” with the desired output file path, while keeping the “${sample}” variable as is at the end of the file path. 

Line 34: module load kallisto/0.46.1
- Change this to load the most recent version of kallisto ONLY IF pseudobam generation is not required. Versions of kallisto following 0.46.1 are not capable of producing pseudobam files. If an updated version of kallisto is to be used then delete lines 72 and 77 to prevent the pseudobam command from executing. 

Line 78: --single -l 180 -s 20 \
- Change the “-l” and “-s” commands based on the fragment lengths and the standard deviation of the fragment lengths of the single end sequencing experiments. 
- THE METHOD TO FIND -L AND -S VALUES ARE IN DEVELOPMENT AND WILL BE ADDED TO THIS REPOSITORY SOON.

To run the sbatch file, the user must use the terminal, navigate to the directory that contains the “Run_Kallisto_pseudobam.sbatch” script, and call upon the script in the following format:

If you are using the default kallisto reference index:
- sbatch Run_Kallisto_pseudobam.sbatch -sample mysamplename 

If you are using a different kallisto reference index:
- sbatch Run_Kallisto_pseudobam.sbatch -sample mysamplename -index indexname

These sbatch calls need to be run for each sample raw data file within the experiment. 

## Part 2: Performing Differential Expression Analysis and Visualization of Results:

The “DESeq2_Analysis.Rmd” script is an R markdown script that performs DESeq2 analysis on the RNA sequenced samples to provide differential expression analysis of genes expressed from different conditions of the experiment. The script first imports the transcript count outputs from the kallisto pseudoalignment step in Part 1, then uses the tx2gene package to convert transcript counts into gene counts. This conversion is required with DESeq2 because the DESeq2 analysis package assumes that the abundance measurements are on the gene-level, not the transcript-level. LFCshrinkage is performed to ensure that the log fold change values are stabilized and not skewed by low count genes in the analysis. To further ensure that the LFC values are stabilized, low count genes are also filtered out before the shrinkage operations. The results of the DESeq2 analysis are then used to generate volcano plots, heat maps, or cluster maps that visualize the results in a clear and readable manner. 

Change the following lines in the “DESeq2_Analysis.Rmd” script based on the needs of the specific RNA-seq experiment. 

Line 27: sample_names <- c('input sample names here')
- Replace the ‘input sample names here’ string with a list of string sample names in which each sample name is within its own set of quotation marks.

Line 29: sample_condition <- c(rep('Condition 1',3),rep('Condition 2',3),rep('Condition 3',3))
- Replace ‘Condition 1’, ‘Condition 2’, and ‘Condition 3’ to reflect the names of the conditions included for the specific RNA-seq experiment. The “rep(‘Condition #,3)” can be added or deleted from this list based on the number of conditions. 

Line 40: analysis_path = "/path/to/kallito/outputs"
- Input the path to the kallisto abundance outputs where it says “path/to/kallisto/outputs”. Make sure that all the outputs are in one directory before running the script. 

Line 44: tx2gene = read.table("/path/to/tx2gene/file/tx2gene.clean.tsv",sep = "\t")
- Input the path to the tx2gene.tsv file where it says “path/to/tx2gene/file/tx2gene.clean.tsv”. 

To create a tx2gene.tsv file, execute the following command in the terminal and within the directory containing the reference fasta sequence used to generate the kallisto index in Part 1:

grep "^>" human_reference.fa | sed -E 's/^>([^ ]+).*gene:([^ ]+).*/\1\t\2/' | awk '{if(NF==1){print $1"\t"$1} else {print $0}}' > tx2gene.tsv

The “human_reference.fa” can be any reference fasta sequence that was also used to generate the kallisto index. This will output a “tx2gene.tsv” file in which there are two columns where the first column has all transcript IDs from the reference, and the second column has the corresponding gene ID of the transcripts from the first column. 

Below is an extra step that may be unnecessary if the first command executes properly. This command will delete any “>” symbols from both columns of the “tx2gene.tsv” file and output it as “tx2gene.clean.tsv”. Any “>” remaining in the file will cause errors within the DESeq2 analysis. 

awk '{gsub(/^>/,"",$1); gsub(/^>/,"",$2); print $1 "\t" $2}' tx2gene.tsv > tx2gene.clean.tsv

Line 190: res_plot$significant = ifelse(!is.na(res_plot$padj) & res_plot$padj < 0.05, "Significant", "Not Significant")  
- Change the 0.05 p-value significance threshold if needed. This value will be the p-value threshold of statistical significance for the volcano plots. 

Call the “volcano_plot” function with the proper input variables for the specific conditions you want visualized. 
- Set the “res” input variable as the specific DESeq2 results object with the desired conditions for the volcano plot. If there are only two conditions, there should only be one result object. If there are more than two conditions, then there should be a list of result objects that can be indexed and called for this input variable. 
- Set the “lfc_val” input variable as the desired threshold value for LFC for a differentially expressed gene to be considered significant. 
- Set the “comp_name” input variable as a string that represents the comparison of conditions for this volcano plot. This string will be included in the title of the plot.

Call the “cluster_map” function with the proper input variables to generate a clustermap with genes of interest on the y-axis, and the sample names on the x-axis. The expression values of the genes are also displayed on the squares of the clustermap.
- Set the “dds” input variable as the specific DESeq2 object with the desired conditions for the clustermap. If there are only two conditions, there should only be one dds object. If there are more than two conditions, then there should be a list of dds objects that can be indexed and called for this input variable. 
- Set the “genes” input variable as a string list of gene symbols that should be included in the cluster map. This can be a collection of gene symbols known to be pluripotency markers. Make sure to define this list of gene symbols as a separate list variable, then setting that list as this input variable. 
- Keep the “gene_names” input variable as “gene_names”. 
- Set the “plot_name” input variable as a string that represents the comparison of conditions for this cluster map. This string will be included in the title of the plot.
  
Call the “int_heatmap” function with the proper input variables to generate a heatmap with genes of interest on the y-axis, and the sample names on the x-axis. This plot is very similar to the cluster map, but either function can be used based on the user’s personal preference.
- The input variables follow the same rules as the input variables for the cluster map. 

