#############################################################
### Get sample names to pass them into a yaml file easily ###
#############################################################

#### wd
setwd("/imppc/labs/eclab/ijarne/0_Recerca/1_NF2_exon11/1_RNA_Seq/2_Gemma_samples_batch2_timepoints/0_RNA_Seq_pipeline")

#### get files (get just R1 and then the pipeline will be able to get R2)
list_of_samples <- list.files(path = "./0_fastq_concatenation/", pattern = "-idt-UMI_1_concat.fastq.gz", full.names = FALSE)
list_of_samples <- gsub(x = list_of_samples, pattern = "-idt-UMI_1_concat.fastq.gz", replacement = "")
writeLines(text = list_of_samples, "./config/samples_list.txt")

