import pandas as pd
from shutil import copyfile
import os, glob, re,sys

##################################
#### Config file information ####
##################################
### Locate the config file
configfile: "./config/config.yaml"
## Access the config file data
samples=config["samples"]
# important paths
base_path = config["base_path"]
# work data
genome_index = config["genome_index"]
gtf_file = config["gft_file"]
adapters = config["adapters"]
# programs
bbduk_sh = config["bbduk_sh"]
# images
fCounts = config["fCounts"]


rule all:
    input:
        ## fastQC htmls and zip of all samples (R1 and R2)
        expand(base_path + "/1_quality_control_after_concatenation/{sample}-idt-UMI_{frr}_concat_fastqc.{extension}",
        sample=config["samples"], frr=["1","2"], extension=["zip","html"],base_path=base_path),
        ## Trimming
        expand(base_path + "/2_trimming/{sample}-idt-UMI_{frr}_trimmed.fastq.gz",
        sample=config["samples"], frr=["1","2"], base_path = base_path),
        ## Post quality trimming
        expand(base_path + "/2_trimming/z_quality_control/{sample}-idt-UMI_{frr}_trimmed_fastqc.{extension}",
        sample=config["samples"], frr=["1","2"], extension=["zip","html"],base_path=base_path),
        ## MultiQC
        base_path + "/1_quality_control_after_concatenation/multiqc_report.html",
        base_path + "/2_trimming/z_quality_control/multiqc_report.html",
        ## Alignment
        # expand(base_path + "/3_alignment/{sample}_Aligned.sortedByCoord.out.bam",
        # sample=config["samples"]),
        ## feature Counts
        base_path + "/4_featureCounts/featureCounts.txt"

        
### Rule FastQC pre trimming       
rule fastQC_first:
    input:
        fastq_1 = base_path + "/0_fastq_concatenation/{sample}-idt-UMI_1_concat.fastq.gz",
        fastq_2 = base_path + "/0_fastq_concatenation/{sample}-idt-UMI_2_concat.fastq.gz"
    output:
        qual_fastq_1_html = base_path + "1_quality_control_after_concatenation/{sample}-idt-UMI_1_concat_fastqc.html",
        qual_fastq_1_zip = base_path + "1_quality_control_after_concatenation/{sample}-idt-UMI_1_concat_fastqc.zip",
        qual_fastq_2_html = base_path + "1_quality_control_after_concatenation/{sample}-idt-UMI_2_concat_fastqc.html",
        qual_fastq_2_zip = base_path + "1_quality_control_after_concatenation/{sample}-idt-UMI_2_concat_fastqc.zip"
    params:
        out_folder = base_path + "/1_quality_control_after_concatenation/"
    threads:
        4
    shell:"""
        fastqc {input.fastq_1} -o {params.out_folder} -t {threads}
        fastqc {input.fastq_2} -o {params.out_folder} -t {threads}
        """
### Rule trimming       
rule bbduk_trimming:
    input:
        fastq_1 = base_path + "/0_fastq_concatenation/{sample}-idt-UMI_1_concat.fastq.gz",
        fastq_2 = base_path + "/0_fastq_concatenation/{sample}-idt-UMI_2_concat.fastq.gz"
    output:
        fastq_1_trimmed = base_path + "/2_trimming/{sample}-idt-UMI_1_trimmed.fastq.gz",
        fastq_2_trimmed = base_path + "/2_trimming/{sample}-idt-UMI_2_trimmed.fastq.gz"
    params:
        adapters = adapters,
        stats = base_path + "/2_trimming/{sample}_trimming_stats.txt",
        refstats = base_path + "/2_trimming/{sample}_trimming_refstats.txt",
        bbduk_sh = bbduk_sh
    threads:
        2
    shell:"""
        {params.bbduk_sh} -Xmx1g in1={input.fastq_1} in2={input.fastq_2} out1={output.fastq_1_trimmed} out2={output.fastq_2_trimmed} ref={params.adapters} ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=20 minlen=23 stats={params.stats} refstats={params.refstats}
        """

### Rule FastQC post trimming       
rule fastQC_post_trimming:
    input:
        fastq_1_trimmed = base_path + "/2_trimming/{sample}-idt-UMI_1_trimmed.fastq.gz",
        fastq_2_trimmed = base_path + "/2_trimming/{sample}-idt-UMI_2_trimmed.fastq.gz"
    output:
        qual_fastq_1_html = base_path + "/2_trimming/z_quality_control/{sample}-idt-UMI_1_trimmed_fastqc.html",
        qual_fastq_1_zip = base_path + "/2_trimming/z_quality_control/{sample}-idt-UMI_1_trimmed_fastqc.zip",
        qual_fastq_2_html = base_path + "/2_trimming/z_quality_control/{sample}-idt-UMI_2_trimmed_fastqc.html",
        qual_fastq_2_zip = base_path + "/2_trimming/z_quality_control/{sample}-idt-UMI_2_trimmed_fastqc.zip"
    params:
        out_folder = base_path + "/2_trimming/z_quality_control/",
        threads_per_sample = 2
    threads:
        4
    shell:"""
        fastqc {input.fastq_1_trimmed} -o {params.out_folder} -t {params.threads_per_sample}
        fastqc {input.fastq_2_trimmed} -o {params.out_folder} -t {params.threads_per_sample}
        """

### MultiQC
rule multiQC_all_samples:
    input:
        fastqc_pre = expand(base_path + "/1_quality_control_after_concatenation/{sample}-idt-UMI_{frr}_concat_fastqc.{extension}",
        sample=config["samples"], frr=["1","2"], extension=["zip","html"],base_path=base_path),
        fastqc_post = expand(base_path + "/2_trimming/z_quality_control/{sample}-idt-UMI_{frr}_trimmed_fastqc.{extension}",
        sample=config["samples"], frr=["1","2"], extension=["zip","html"],base_path=base_path)
    output:
        html_pre = base_path + "/1_quality_control_after_concatenation/multiqc_report.html",
        html_post = base_path + "/2_trimming/z_quality_control/multiqc_report.html"
    params:
        out_pre = base_path + "/1_quality_control_after_concatenation/",
        out_post = base_path + "/2_trimming/z_quality_control/"
    threads:
        2
    shell:"""
        multiqc {params.out_pre} --outdir {params.out_pre} 
        multiqc {params.out_post} --outdir {params.out_post}
        """

### Rule STAR alignment
# rule star_alignment:
#     input:
#         fastq_1_trimmed = base_path + "/2_trimming/{sample}-idt-UMI_1_trimmed.fastq.gz",
#         fastq_2_trimmed = base_path + "/2_trimming/{sample}-idt-UMI_2_trimmed.fastq.gz"
#     output:
#         aligned = base_path + "/3_alignment/{sample}_Aligned.sortedByCoord.out.bam"
#     params:
#         genome_index = genome_index,
#         threads = 12,
#         out_prefix = base_path + "/3_alignment/{sample}_"
#     threads:
#         12
#     shell:"""
#         STAR --runThreadN {params.threads} --genomeDir {params.genome_index} \
#         --readFilesIn {input.fastq_1_trimmed} {input.fastq_2_trimmed} --readFilesCommand zcat \
#         --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN {params.threads} \
#         --outFileNamePrefix {params.out_prefix}
#         """

### Rule featureCounts
rule counting:
    input:
        bams_to_count = expand(base_path + "/3_alignment/{sample}_Aligned.sortedByCoord.out.bam",
        sample=config["samples"])
    output:
        fc_table = base_path + "/4_featureCounts/featureCounts.txt"
    params:
        gtf_file = gtf_file,
    threads:
        10
    container:
        fCounts
    shell:"""
        featureCounts -T 10 -t exon -g gene_id \
        --primary -p -s 1  -C -B -a {params.gtf_file} \
        -o {output.fc_table} \
        {input.bams_to_count}
        """

# STAR --runThreadN 12 --genomeDir "/imppc/labs/eclab/ijarne/0_Recerca/1_NF2_exon11/1_RNA_Seq/0_Nuria_work/0_RNA_seq_pipeline/2_alignment/0_index_preparation/indexed_genome_STAR" \
# --readFilesIn "${pair}_1_trimmed.fastq.gz" "${pair}_2_trimmed.fastq.gz" --readFilesCommand zcat \
# --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 12\
# --outFileNamePrefix "${CURRENT_DIRECTORY}/2_alignment/${pair}"



### Rule MultiQC       
# rule multiqc_before_trimming:
#     input:
#         expand(base_path + "/0_quality_metrics/0_before_trimming/{sample}_L001_{frr}_001_fastqc.{extension}",
#         sample=config["samples"], frr=["R1","R2"], extension=["zip","html"],base_path=base_path)
#     output:
#         base_path + "/0_quality_metrics/0_before_trimming/report.html"
#     threads:
#         4
#     params:
#         out_name = "report.html",
#         out_route = base_path + "/0_quality_metrics/0_before_trimming/"
#     conda:
#         anac_envs + config["Anaconda_envs"]["multiqc"]
#     benchmark:
#         base_path + "/Run/benchmarks/beforetrimming_multiQC.bmk"
#     log:
#         base_path + "Run/logs/beforetrimming_multiQC.log"
#     shell:"""
#         multiqc --filename {params.out_name} --outdir {params.out_route} {params.out_route} 2>> {log}
#         """

