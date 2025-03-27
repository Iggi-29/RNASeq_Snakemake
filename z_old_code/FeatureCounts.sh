#!/bin/bash
# request Bourne shell as shell for job
#$ -S /bin/bash
# Execute from the current workig dir
#$ -cwd
# Name for the script in the queuing system
#$ -N batch_NF2_fCounts_v2_No_PAIRS
# In order to load the environment variables and your path
# You can either use this or do a : source /etc/profile
#$ -V
# You can redirect the error output to a specific file
#$ -e batch_NF2_fCounts_v2_No_PAIRS.err
# You can redirect the output to a specific file
#$ -o batch_NF2_fCounts_v2_No_PAIRS.log
#$ -pe smp 12
#$ -q d10imppcv3 # los nodos nuevos!
#$ -l h_vmem=9G
# Avoid node 16 (low memory at the moment)
# -l h=!sge-exec-16

echo "#############################"
echo "######  featureCouns   ######"
echo "#############################"

source /imppc/labs/eclab/ijarne/miniconda3/etc/profile.d/conda.sh
conda activate base

### Save the names of some folders that can be useful for later
CURRENT_DIRECTORY=$(pwd)

singularity exec -e -B /imppc/labs/eclab/ "${CURRENT_DIRECTORY}/featurecounts_latest.sif" featureCounts -T 10 -t exon -g gene_id --primary -p -s 1  -C -B -a "${CURRENT_DIRECTORY}/3_featureCounts/Homo_sapiens.GRCh38.113.gtf.gz" -o ./3_featureCounts/featureCounts_C_B.txt  ./3_featureCounts/sorted_bams/*.sortedByCoord.out.bam

echo "Finished !"