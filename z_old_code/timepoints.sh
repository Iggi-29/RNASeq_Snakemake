#!/bin/bash
# request Bourne shell as shell for job
#$ -S /bin/bash
# Execute from the current workig dir
#$ -cwd
# Name for the script in the queuing system
#$ -N batch_time
# In order to load the environment variables and your path
# You can either use this or do a : source /etc/profile
#$ -V
# You can redirect the error output to a specific file
#$ -e batch_time.err
# You can redirect the output to a specific file
#$ -o batch_time.log
#$ -pe smp 16
#$ -q d10imppcv3 # los nodos nuevos!
#$ -l h_vmem=9G
# Avoid node 16 (low memory at the moment)
# -l h=!sge-exec-16


echo "#############################"
echo "#########  fatsqc   #########"
echo "#############################"
### Save the names of some folders that can be useful for later
RAW_DATA_1="/imppc/labs/eclab/ijarne/0_Recerca/1_NF2_exon11/1_RNA_Seq/2_Gemma_samples_batch2_timepoints/0_RNA_Seq_pipeline/0_fastq_concatenation" # Place of the Raw fastq files
CURRENT_DIRECTORY=$(pwd) # Place where we want to save the files
cd $RAW_DATA_1

echo "The data is stored here: $RAW_DATA_1"
echo "The report will be at $CURRENT_DIRECTORY"

counter_1=0
### Get the info from the name of each file
for file in *.fastq.gz; do
   counter_1=$((counter_1 +1))
   echo "Processing file #$counter_1: $file"
   fastqc "${file}" -o "${CURRENT_DIRECTORY}/1_quality_control_after_concatenation" -t 18
   echo "Finished with pair #$counter_1: $pair"
done

echo "Fastq Finished"

echo "###########################"
echo "#########  bbduk  #########"
echo "###########################"
### Save the names of some folders that can be useful for later
RAW_DATA_2="/imppc/labs/eclab/ijarne/0_Recerca/1_NF2_exon11/1_RNA_Seq/2_Gemma_samples_batch2_timepoints/0_RNA_Seq_pipeline/0_fastq_concatenation" # Place of the Raw fastq files
cd $RAW_DATA_2

## Save bbduk file locations
bbduk="/imppc/labs/eclab/ijarne/miniconda3/opt/bbmap-39.18-0/bbduk.sh"
adapters="/imppc/labs/eclab/ijarne/miniconda3/opt/bbmap-39.18-0/resources/adapters.fa"

counter_2=0
### Get the info from the name of each file
for pair in `ls *_1_concat.fastq.gz | sed 's/_1_concat.fastq.gz//'`; do
    counter_2=$((counter_2 +1))
    echo "Processing pair #$counter_2: $pair"
    $bbduk -Xmx1g in1="${pair}_1_concat.fastq.gz" in2="${pair}_2_concat.fastq.gz" out1="${CURRENT_DIRECTORY}/2_trimming/${pair}_1_trimmed.fastq.gz" out2="${CURRENT_DIRECTORY}/2_trimming/${pair}_2_trimmed.fastq.gz" ref="${adapters}" ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=20 minlen=23 stats="${CURRENT_DIRECTORY}/2_trimming/${pair}_trimming_stats.txt" refstats="${CURRENT_DIRECTORY}/2_trimming/${pair}_trimming_refstats.txt"
    echo "Finished with pair #$counter_2: $pair"
done  

echo "###############################"
echo "#########  fatsqc_2   #########"
echo "###############################"
### Save the names of some folders that can be useful for later
RAW_DATA_3="/imppc/labs/eclab/ijarne/0_Recerca/1_NF2_exon11/1_RNA_Seq/2_Gemma_samples_batch2_timepoints/0_RNA_Seq_pipeline/2_trimming" # Place of the Raw fastq files
cd $RAW_DATA_3

counter_3=0
for file in *.fastq.gz; do
    counter_3=$((counter_3 +1))
    echo "Processing file #$counter_3: $file"
    fastqc "${file}" -o "${CURRENT_DIRECTORY}/2_trimming/z_quality_control/" -t 14
    echo "Finished with file #$counter_3: $file"
done

echo "###########################"
echo "#########  STAR   #########"
echo "##########################"

### Save the names of some folders that can be useful for later
RAW_DATA_4="/imppc/labs/eclab/ijarne/0_Recerca/1_NF2_exon11/1_RNA_Seq/2_Gemma_samples_batch2_timepoints/0_RNA_Seq_pipeline/2_trimming" # Place of the Raw fastq files
cd $RAW_DATA_4

counter_4=0
for pair in `ls *_1_trimmed.fastq.gz | sed 's/_1_trimmed.fastq.gz//'`; do
    counter_4=$((counter_4 +1))
    echo "Processing pair #$counter_4: $pair"
    STAR --runThreadN 12 --genomeDir "/imppc/labs/eclab/ijarne/0_Recerca/1_NF2_exon11/1_RNA_Seq/0_Nuria_work/0_RNA_seq_pipeline/2_alignment/0_index_preparation/indexed_genome_STAR" \
    --readFilesIn "${pair}_1_trimmed.fastq.gz" "${pair}_2_trimmed.fastq.gz" --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 12\
    --outFileNamePrefix "${CURRENT_DIRECTORY}/2_alignment/${pair}"
    echo "Finished with pair #$counter_4: $pair"
done 

