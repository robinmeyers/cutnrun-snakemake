# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.

import re
import gzip
import glob
import csv
import pandas as pd
from snakemake.utils import validate, min_version
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


##### load config and sample sheets #####


configfile: "config.yaml"
# report: "report/workflow.rst"

validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_csv(config['samplesheet']).set_index("Name", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")


SAMPLES = dict(zip(samples["Name"], samples["ID"]))

THREADS = config['threads']
FASTQ_DIR = config['fastq_dir']



# Allow users to fix the underlying OS via singularity.
singularity: "docker://continuumio/miniconda3"


# include: "rules/other.smk"

wildcard_constraints:
    directory=".+\/",
    sample="[^\/\.]+"


localrules: all


rule all:
    input:
        bigwigs = expand("outs/signal/{sample}.bigwig", sample=SAMPLES.keys()),
        peaks = expand("outs/peaks/{sample}.stringent.bed", sample=SAMPLES.keys())
    run:
        print("workflow complete!")


def get_paired_fqs(wildcards):
    sample_id = SAMPLES[wildcards.sample]
    r1 = glob.glob(os.path.join(FASTQ_DIR, "**", sample_id + "_*R1_*.fastq.gz"),
        recursive=True)
    r2 = glob.glob(os.path.join(FASTQ_DIR, "**", sample_id + "_*R2_*.fastq.gz"), 
        recursive=True)
    if len(r1) == 0:
        raise ValueError(sample_id + " has no matching input fastq file")
    if len(r1) != len(r2):
        raise ValueError(sample_id + " has different numbers of R1 and R2 fastq files")
    return {"r1": sorted(r1), "r2": sorted(r2)}

rule run_fastqc:
    input: unpack(get_paired_fqs)
    output: "outs/fastqc/{sample}/.done"
    shell:
        "fastqc -o outs/fastqc/{wildcards.sample} {input.r1} {input.r2} && touch {output}"


rule merge_fastqs:
    input: unpack(get_paired_fqs)
    output:
        r1 = temp("outs/merge/{sample}_R1.fastq.gz"),
        r2 = temp("outs/merge/{sample}_R2.fastq.gz")
    shell:
        "cat {input.r1} > {output.r1}; cat {input.r2} > {output.r2}"


rule trim_adaptors:
    input: 
        r1 = "outs/merge/{sample}_R1.fastq.gz",
        r2 = "outs/merge/{sample}_R2.fastq.gz"
    output:
        r1 = "outs/trim/{sample}_R1.fastq.gz",
        r2 = "outs/trim/{sample}_R2.fastq.gz"
        # qc = "outs/trim/{sample}_qc.txt"
    threads: THREADS
    log: "outs/trim/{sample}.log"
    shell:
        "cutadapt -a {config[adaptor_5p]} -A {config[adaptor_3p]} -m {config[min_len]} --cores {threads} "
        "-o {output.r1} -p {output.r2} {input.r1} {input.r2} > {log}"

rule align:
    input:
        r1 = "outs/trim/{sample}_R1.fastq.gz",
        r2 = "outs/trim/{sample}_R2.fastq.gz"
    output:
        sam = "outs/align/{sample}.sam",
        log = "outs/align/{sample}.bowtie2.log"
    threads: THREADS
    shell:
        "bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 "
        "-p {threads} -x {config[reference]} "
        "-1 {input.r1} -2 {input.r2} -S {output.sam} &> {output.log}"


rule sort_filter_bam:
    input:
        "outs/align/{sample}.bam"
    output:
        "outs/align/{sample}.cleaned.bam"
    threads: THREADS
    shell:
        "samtools fixmate -m -@ {threads} {input} - | "
        "samtools sort -@ {threads}  - | "
        "samtools markdup - - | "
        "samtools view -b -f 0x3 -F 0x400 - | "
        "samtools sort -n - > {output}"
        

rule bam_to_bed:
    input: "outs/align/{sample}.cleaned.bam"
    output: "outs/align/{sample}.bed"
    shell:
        "bedtools bamtobed -bedpe -i {input} > {output}; "
  

rule clean_bed:
    input: "outs/align/{sample}.bed"
    output: "outs/align/{sample}.cleaned.bed"
    shell:
        "awk '$1==$4 && $6-$2 < 1000 {{print $0}}' {input} | "
        "cut -f 1,2,6 | sort -k1,1 -k2,2n -k3,3n > {output} "

rule bed_to_bedgraph:
    input: "outs/align/{sample}.cleaned.bed"
    output: "outs/signal/{sample}.bedgraph"
    shell:
        "bedtools genomecov -bg -i {input} -g {config[chrom_sizes]} > {output}"

rule bedgraph_to_bigwig:
    input: "outs/signal/{sample}.bedgraph"
    output: "outs/signal/{sample}.bigwig"
    shell:
        "bedGraphToBigWig {input} {config[chrom_sizes]} {output}"


rule call_peaks:
    input: "outs/signal/{sample}.bedgraph"
    output: "outs/peaks/{sample}.stringent.bed"
    shell:
        "SEACR_1.3.sh {input} 0.1 norm stringent outs/peaks/{wildcards.sample}"



# rule mapping_stats:
#     input:
#     output:
#     shell:




rule sam_to_bam:
    input:
        "{directory}{sample}.sam"
    output:
        "{directory}{sample}.bam"
    threads: THREADS
    shell:
        "samtools view -b --threads {threads} -o {output} {input}"




    
rule samtools_index:
    input:
        "{directory}{sample}.bam"
    output:
        "{directory}{sample}.bam.bai"
    threads: THREADS
    shell:
        "samtools index -@ {threads} {input}"
