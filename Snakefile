# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.

import re
import gzip
import glob
import csv
import numpy as np
import pandas as pd
from snakemake.utils import validate, min_version
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


##### load config and sample sheets #####


configfile: "config.yaml"
# report: "report/workflow.rst"

validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_csv(config['samplesheet']).set_index("sample", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")


SAMPLES = dict(zip(samples["sample"], samples["fastq"]))


if "control" in samples.columns:
    samples_with_controls = samples[samples["control"].notnull()]
    CONTROLS = dict(zip(samples_with_controls["sample"], samples_with_controls["control"]))


if "condition" in samples.columns:
    CONDITIONS = {}
    for condition in np.unique(samples["condition"]):
        CONDITIONS[condition] = np.unique(samples[samples["condition"] == condition]["sample"])
    if "control" in samples.columns:
        samples_with_controls = samples[samples["control"].notnull()]
        CONDITION_CONTROLS = {}
        for condition in np.unique(samples_with_controls["condition"]):
            CONDITION_CONTROLS[condition] = np.unique(samples[samples["condition"] == condition]["control"])




THREADS = config['threads']
FASTQ_DIR = config['fastq_dir']





# Allow users to fix the underlying OS via singularity.
singularity: "docker://continuumio/miniconda3"


# include: "rules/other.smk"

wildcard_constraints:
    directory=".+\/",
    filename="[^\/]+",
    sample="[^\/\.]+",
    condition="[^\/\.]+",
    dir_type="(samples|conditions)",
    stringency="(relaxed|stringent)"


localrules: all


def get_target_files(wildcards):
    targets = []
    targets = targets + expand("outs/samples/fastqc/{sample}/.done", sample=SAMPLES.keys())

    targets = targets + expand("outs/samples/signal/{sample}.bigwig", sample=SAMPLES.keys())
    targets = targets + expand("outs/samples/peaks/{sample}.{stringency}.bed", sample=SAMPLES.keys(), stringency = ['relaxed', 'stringent'])
    if "control" in samples.columns:
        targets = targets + expand("outs/samples/peaks/{sample}.vs-ctrl.{stringency}.bed", sample=CONTROLS.keys(), stringency = ['relaxed', 'stringent'])
    if "condition" in samples.columns:
        targets = targets + expand("outs/conditions/signal/{condition}.bigwig", condition=CONDITIONS.keys())
        targets = targets + expand("outs/conditions/peaks/{condition}.{stringency}.bed", condition=CONDITIONS.keys(), stringency = ['relaxed', 'stringent'])
        if "control" in samples.columns:
            targets = targets + expand("outs/conditions/peaks/{condition}.vs-ctrl.{stringency}.bed", condition=CONDITION_CONTROLS.keys(), stringency = ['relaxed', 'stringent'])
    if config['reference_spikein']:
        targets = targets + expand("outs/samples/align-spikein/{sample}.spikein.bam", sample=SAMPLES.keys())
    return targets

rule all:
    input: get_target_files
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
    output: "outs/samples/fastqc/{sample}/.done"
    threads: 1
    shell:
        "fastqc -o outs/samples/fastqc/{wildcards.sample} {input.r1} {input.r2} && touch {output}"


rule merge_fastqs:
    input: unpack(get_paired_fqs)
    output:
        r1 = temp("outs/samples/merge/{sample}_R1.fastq.gz"),
        r2 = temp("outs/samples/merge/{sample}_R2.fastq.gz")
    threads: 1
    shell:
        "cat {input.r1} > {output.r1}; cat {input.r2} > {output.r2}"


rule trim_adaptors:
    input: 
        r1 = "outs/samples/merge/{sample}_R1.fastq.gz",
        r2 = "outs/samples/merge/{sample}_R2.fastq.gz"
    output:
        r1 = "outs/samples/trim/{sample}_R1.fastq.gz",
        r2 = "outs/samples/trim/{sample}_R2.fastq.gz"
        # qc = "outs/trim/{sample}_qc.txt"
    threads: THREADS
    log: "outs/samples/trim/{sample}.log"
    shell:
        "cutadapt -a {config[adaptor_5p]} -A {config[adaptor_3p]} -m {config[min_len]} --cores {threads} "
        "-o {output.r1} -p {output.r2} {input.r1} {input.r2} > {log}"

rule align:
    input:
        r1 = "outs/samples/trim/{sample}_R1.fastq.gz",
        r2 = "outs/samples/trim/{sample}_R2.fastq.gz"
    output:
        sam = "outs/samples/align/{sample}.sam",
        log = "outs/samples/align/{sample}.bowtie2.log"
    threads: THREADS
    shell:
        "bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 "
        "-p {threads} -x {config[reference]} "
        "-1 {input.r1} -2 {input.r2} -S {output.sam} &> {output.log}"

rule align_spikein:
    input:
        r1 = "outs/samples/trim/{sample}_R1.fastq.gz",
        r2 = "outs/samples/trim/{sample}_R2.fastq.gz"
    output:
        sam = "outs/samples/align-spikein/{sample}.spikein.sam",
        log = "outs/samples/align-spikein/{sample}.spikein.bowtie2.log"
    threads: THREADS
    shell:
        "bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 --no-unal "
        "-p {threads} -x {config[reference_spikein]} "
        "-1 {input.r1} -2 {input.r2} -S {output.sam} &> {output.log}"


rule sort_filter_bam:
    input:
        "outs/samples/align/{sample}.bam"
    output:
        "outs/samples/align/{sample}.cleaned.bam"
    threads: THREADS
    shell:
        "samtools fixmate -m -@ {threads} {input} - | "
        "samtools sort -@ {threads} -T outs/samples/align/{wildcards.sample} - | "
        "samtools markdup - - | "
        "samtools view -b -f 0x3 -F 0x400 - | "
        "samtools sort -n - > {output}"
        

rule bam_to_bed:
    input: "outs/samples/align/{sample}.cleaned.bam"
    output: "outs/samples/align/{sample}.bed"
    threads: 1
    shell:
        "bedtools bamtobed -bedpe -i {input} > {output}; "
  

rule clean_bed:
    input: "outs/samples/align/{sample}.bed"
    output: "outs/samples/align/{sample}.cleaned.bed"
    threads: 1
    shell:
        "awk '$1==$4 && $6-$2 < 1000 {{print $0}}' {input} | "
        "cut -f 1,2,6 | sort -k1,1 -k2,2n -k3,3n > {output} "

rule bed_to_bedgraph:
    input: "outs/{dir_type}/align/{sample}.cleaned.bed"
    output: "outs/{dir_type}/signal/{sample}.bedgraph"
    threads: 1
    shell:
        "bedtools genomecov -bg -i {input} -g {config[chrom_sizes]} > {output}"

rule bedgraph_to_bigwig:
    input: "outs/{dir_type}/signal/{sample}.bedgraph"
    output: "outs/{dir_type}/signal/{sample}.bigwig"
    threads: 1
    shell:
        "bedGraphToBigWig {input} {config[chrom_sizes]} {output}"


rule call_seacr_peaks:
    input: "outs/{dir_type}/signal/{sample}.bedgraph"
    output: "outs/{dir_type}/peaks/{sample}.{stringency}.bed"
    threads: 1
    shell:
        "SEACR_1.3.sh {input} 0.01 non {wildcards.stringency} "
        "outs/{wildcards.dir_type}/peaks/{wildcards.sample}.{wildcards.stringency}; "
        "mv outs/{wildcards.dir_type}/peaks/{wildcards.sample}.{wildcards.stringency}.{wildcards.stringency}.bed "
        "outs/{wildcards.dir_type}/peaks/{wildcards.sample}.{wildcards.stringency}.bed"


def get_expt_and_ctrl_bedgraphs(wildcards):
    if wildcards.dir_type == "samples":
        control_sample = CONTROLS[wildcards.sample]
    if wildcards.dir_type == "conditions":
        control_sample = wildcards.sample + "_CONTROL"
    expt_bg = os.path.join("outs", wildcards.dir_type, "signal", wildcards.sample + ".bedgraph")
    ctrl_bg = os.path.join("outs", wildcards.dir_type, "signal", control_sample + ".bedgraph")
    return {'expt' : expt_bg, 'ctrl' : ctrl_bg}

rule call_seacr_peaks_vs_control:
    input: unpack(get_expt_and_ctrl_bedgraphs)
    output: "outs/{dir_type}/peaks/{sample}.vs-ctrl.{stringency}.bed"
    threads: 1
    shell:
        "SEACR_1.3.sh {input.expt} {input.ctrl} norm {wildcards.stringency} "
        "outs/{wildcards.dir_type}/peaks/{wildcards.sample}.vs-ctrl.{wildcards.stringency}; "
        "mv outs/{wildcards.dir_type}/peaks/{wildcards.sample}.vs-ctrl.{wildcards.stringency}.{wildcards.stringency}.bed "
        "outs/{wildcards.dir_type}/peaks/{wildcards.sample}.vs-ctrl.{wildcards.stringency}.bed"


def get_samples_per_condition(wildcards):
    return [os.path.join("outs/samples/align/", s + ".cleaned.bed") for s in CONDITIONS[wildcards.condition]]

rule merge_samples:
    input: get_samples_per_condition
    output: "outs/conditions/align/{condition}.cleaned.bed"
    threads: 1
    shell:
        "cat {input} | sort -k1,1 -k2,2n -k3,3n > {output}"

def get_controls_per_condition(wildcards):
    return [os.path.join("outs/samples/align/", s + ".cleaned.bed") for s in CONDITION_CONTROLS[wildcards.condition]]

rule merge_controls:
    input: get_controls_per_condition
    output: "outs/conditions/align/{condition}_CONTROL.cleaned.bed"
    threads: 1
    shell:
        "cat {input} | sort -k1,1 -k2,2n -k3,3n > {output}"



# rule mapping_stats:
#     input:
#     output:
#     shell:




rule sam_to_bam:
    input:
        "{directory}{filename}.sam"
    output:
        "{directory}{filename}.bam"
    threads: THREADS
    shell:
        "samtools view -b --threads {threads} -o {output} {input}"




    
rule samtools_index:
    input:
        "{directory}{filename}.bam"
    output:
        "{directory}{filename}.bam.bai"
    threads: THREADS
    shell:
        "samtools index -@ {threads} {input}"
