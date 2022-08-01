# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.

import re
# import gzip
import glob
# import csv
import numpy as np
import pandas as pd
from snakemake.utils import validate, min_version
# from Bio import SeqIO
# from Bio.SeqRecord import SeqRecord
# from Bio.Seq import Seq

# report: "report/workflow.rst"

##### load config and sample sheets #####
configfile: "config.yaml"
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
OUTPUT_DIR = os.path.join(config['output_dir'], '') if config['output_dir'] else "outs/"


def prefix_out(filename):
    return(os.path.join(OUTPUT_DIR, filename))

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
    targets = targets + expand("samples/fastqc/{sample}/.done", sample=SAMPLES.keys())

    targets = targets + expand("samples/signal/{sample}.scaled.bigwig", sample=SAMPLES.keys())

    
    if "control" in samples.columns:
        targets = targets + expand("samples/peaks/{sample}.vs-ctrl.{stringency}.filtered.bed", sample=CONTROLS.keys(), stringency = ['relaxed', 'stringent'])
        targets = targets + expand("samples/peaks_macs2/{sample}.macs2_q0.1_peaks.filtered.narrowPeak", sample=CONTROLS.keys())
    else:
        targets = targets + expand("samples/peaks/{sample}.{stringency}.filtered.bed", sample=SAMPLES.keys(), stringency = ['relaxed', 'stringent'])
    if "condition" in samples.columns:
        targets = targets + expand("conditions/signal/{condition}.scaled.bigwig", condition=CONDITIONS.keys())

        if "control" in samples.columns:
            targets = targets + expand("conditions/peaks/{condition}.vs-ctrl.{stringency}.filtered.bed", condition=CONDITION_CONTROLS.keys(), stringency = ['relaxed', 'stringent'])
            targets = targets + expand("conditions/peaks_macs2/{sample}.macs2_q0.1_peaks.filtered.narrowPeak", sample=CONDITION_CONTROLS.keys())

        else:
            targets = targets + expand("conditions/peaks/{condition}.{stringency}.filtered.bed", condition=CONDITIONS.keys(), stringency = ['relaxed', 'stringent'])

    # if config['reference_spikein']:
    #     targets = targets + expand("samples/align-spikein/{sample}.spikein.bam.seqdepth", sample=SAMPLES.keys())

    targets = targets + expand("samples/align/{sample}.cleaned.fragmentsize.txt", sample = SAMPLES.keys())
    targets = targets + ["samples/alignment_summary.csv"]

    targets = [os.path.join(OUTPUT_DIR, t) for t in targets]
    return targets

rule all:
    input: get_target_files
    run:
        print("workflow complete!")


def get_paired_fqs(wildcards):
    if (wildcards.sample not in SAMPLES.keys()):
        return {'r1' : [], 'r2': []}
        # raise ValueError(wildcards.sample + " is not a sample in the samplesheet")
    sample_id = SAMPLES[wildcards.sample]

    r1 = list(filter(re.compile(sample_id + "(_S[0-9]+)?(_L[0-9]+)?_R1(_001)?.fastq.gz$").search, 
        glob.glob(os.path.join(FASTQ_DIR, sample_id + "*"), recursive=True)))
    r2 = list(filter(re.compile(sample_id + "(_S[0-9]+)?(_L[0-9]+)?_R2(_001)?.fastq.gz$").search, 
        glob.glob(os.path.join(FASTQ_DIR, sample_id + "*"), recursive=True)))

    # r1 = glob.glob(os.path.join(FASTQ_DIR, "**", sample_id + "_*R1*.fastq.gz"),
    #     recursive=True)
    # r2 = glob.glob(os.path.join(FASTQ_DIR, "**", sample_id + "_*R2*.fastq.gz"), 
    #     recursive=True)
    if len(r1) == 0:
        raise ValueError(sample_id + " has no matching input fastq file")
    if len(r1) != len(r2):
        raise ValueError(sample_id + " has different numbers of R1 and R2 fastq files")
    return {"r1": sorted(r1), "r2": sorted(r2)}


rule run_fastqc:
    input: unpack(get_paired_fqs)
    output: OUTPUT_DIR + "samples/fastqc/{sample}/.done"
    params:
        out = OUTPUT_DIR + "samples/fastqc/{sample}"
    log: OUTPUT_DIR + "samples/fastqc/{sample}/{sample}.log"
    threads: 1
    shell:
        "fastqc -o {params.out} {input.r1} {input.r2} && touch {output}"


rule merge_fastqs:
    input: unpack(get_paired_fqs)
    output:
        r1 = temp(OUTPUT_DIR + "samples/merge/{sample}_R1.fastq.gz"),
        r2 = temp(OUTPUT_DIR + "samples/merge/{sample}_R2.fastq.gz")
    threads: 1
    shell:
        "cat {input.r1} > {output.r1}; cat {input.r2} > {output.r2}"


rule trim_adaptors:
    input: 
        r1 = OUTPUT_DIR + "samples/merge/{sample}_R1.fastq.gz",
        r2 = OUTPUT_DIR + "samples/merge/{sample}_R2.fastq.gz"
    output:
        r1 = OUTPUT_DIR + "samples/trim/{sample}_R1.fastq.gz",
        r2 = OUTPUT_DIR + "samples/trim/{sample}_R2.fastq.gz"
    threads: THREADS
    log: OUTPUT_DIR + "samples/trim/{sample}.log"
    shell:
        "cutadapt -a {config[adaptor_5p]} -A {config[adaptor_3p]} -m {config[min_len]} --cores {threads} "
        "-o {output.r1} -p {output.r2} {input.r1} {input.r2} > {log}"

rule align:
    input:
        r1 = OUTPUT_DIR + "samples/trim/{sample}_R1.fastq.gz",
        r2 = OUTPUT_DIR + "samples/trim/{sample}_R2.fastq.gz"
    output:
        sam = OUTPUT_DIR + "samples/align/{sample}.sam"
    log: OUTPUT_DIR + "samples/align/{sample}.bowtie2.log"
    threads: THREADS
    resources:
        mem_mb = 32000
    shell:
        "bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 "
        "-p {threads} -x {config[reference]} "
        "-1 {input.r1} -2 {input.r2} -S {output.sam} &> {log}"

rule align_spikein:
    input:
        r1 = OUTPUT_DIR + "samples/trim/{sample}_R1.fastq.gz",
        r2 = OUTPUT_DIR + "samples/trim/{sample}_R2.fastq.gz"
    output:
        sam = OUTPUT_DIR + "samples/align-spikein/{sample}.spikein.sam"
    log: OUTPUT_DIR + "samples/align-spikein/{sample}.spikein.bowtie2.log"
    threads: THREADS
    resources:
        mem_mb = 8000
    shell:
        "bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 --no-unal "
        "-p {threads} -x {config[reference_spikein]} "
        "-1 {input.r1} -2 {input.r2} -S {output.sam} &> {log}"


rule seqdepth:
    input: OUTPUT_DIR + "{dir_type}/align/{sample}.cleaned.bam"
    output: OUTPUT_DIR + "{dir_type}/align/{sample}.cleaned.bam.seqdepth"
    shell:
        "samtools view -f 0x42 {input} | wc -l | tr -d ' '> {output}"

rule spikein_seqdepth:
    input: OUTPUT_DIR + "{dir_type}/align-spikein/{sample}.spikein.cleaned.bam"
    output: OUTPUT_DIR + "{dir_type}/align-spikein/{sample}.spikein.cleaned.bam.seqdepth"
    shell:
        "samtools view -f 0x42 {input} | wc -l | tr -d ' '> {output}"

rule fragment_size:
    input: OUTPUT_DIR + "{dir_type}/align/{sample}.cleaned.bam"
    output: OUTPUT_DIR + "{dir_type}/align/{sample}.cleaned.fragmentsize.txt"
    shell:
        "samtools view -f 0x42 {input} | awk -F'\\t' 'function abs(x) {{return ((x < 0.0) ? -x : x)}} {{print abs($9)}}' | "
        "sort -n | uniq -c | awk -v OFS='\\t' '{{print $2, $1}}' > {output}"

rule alignment_summary:
    input:
        lambda wildcards: expand(OUTPUT_DIR + "samples/align/{sample}.cleaned.bam.seqdepth", sample=SAMPLES.keys()),
        lambda wildcards: expand(OUTPUT_DIR + "samples/align-spikein/{sample}.spikein.cleaned.bam.seqdepth", sample=SAMPLES.keys()) if config['reference_spikein'] else []
    output:
        alignment_summary = OUTPUT_DIR + "samples/alignment_summary.csv"
    log: OUTPUT_DIR + "samples/alignment_summary.log"
    script: "scripts/alignment_summary.R"


rule sort_filter_bam:
    input:
        OUTPUT_DIR + "samples/align/{sample}.bam"
    output:
        OUTPUT_DIR + "samples/align/{sample}.cleaned.bam"
    params:
        tmp = OUTPUT_DIR + "samples/align/{sample}"
    threads: THREADS
    resources:
        mem_mb = 32000
    shell:
        "samtools fixmate -m -@ {threads} {input} - | "
        "samtools sort -@ {threads} -T {params.tmp} - | "
        "samtools markdup - - | "
        "samtools view -b -f 0x3 -F 0x400 - | "
        "samtools sort -n - > {output}"
        

rule sort_filter_spikein_bam:
    input:
        OUTPUT_DIR + "samples/align-spikein/{sample}.spikein.bam"
    output:
        OUTPUT_DIR + "samples/align-spikein/{sample}.spikein.cleaned.bam"
    params:
        tmp = OUTPUT_DIR + "samples/align-spikein/{sample}"
    threads: THREADS
    resources:
        mem_mb = 8000
    shell:
        "samtools fixmate -m -@ {threads} {input} - | "
        "samtools sort -@ {threads} -T {params.tmp} - | "
        "samtools markdup - - | "
        "samtools view -b -f 0x3 -F 0x400 - | "
        "samtools sort -n - > {output}"

rule bam_to_bed:
    input: OUTPUT_DIR + "{dir_type}/align/{sample}.cleaned.bam"
    output: OUTPUT_DIR + "{dir_type}/align/{sample}.bed"
    threads: 1
    shell:
        "bedtools bamtobed -bedpe -i {input} > {output}; "
  

rule clean_bed:
    input: OUTPUT_DIR + "{dir_type}/align/{sample}.bed"
    output: OUTPUT_DIR + "{dir_type}/align/{sample}.cleaned.bed"
    threads: 1
    shell:
        "awk '$1==$4 && $6-$2 < 1000 {{print $0}}' {input} | "
        "cut -f 1,2,6,7,8,9 | sort -k1,1 -k2,2n -k3,3n > {output} "

rule bed_to_bedgraph:
    input: OUTPUT_DIR + "{dir_type}/align/{sample}.cleaned.bed"
    output: OUTPUT_DIR + "{dir_type}/signal/{sample}.bedgraph"
    threads: 1
    resources:
        mem_mb = 4000
    shell:
        "bedtools genomecov -bg -i {input} -g {config[chrom_sizes]} > {output}"


# def get_scale_factor(wildcards, input):
#     with open(input.seqdepth) as f:
#         seqdepth = int(f.readline())
#     if config['reference_spikein']:
#         scale_factor = 10000 / seqdepth
#     else:
#         scale_factor = 1000000 / seqdepth
#     return scale_factor

rule bed_to_scaled_bedgraph:
    input:
        bed = OUTPUT_DIR + "{dir_type}/align/{sample}.cleaned.bed",
        seqdepth = lambda wildcards: os.path.join(OUTPUT_DIR, wildcards.dir_type, "align-spikein", wildcards.sample + ".spikein.cleaned.bam.seqdepth") if config['reference_spikein'] else os.path.join(OUTPUT_DIR, wildcards.dir_type, "align", wildcards.sample + ".cleaned.bam.seqdepth")
    output: OUTPUT_DIR + "{dir_type}/signal/{sample}.scaled.bedgraph"
    threads: 1
    params: scale_constant = lambda wildcards: 10000 if config['reference_spikein'] else 10000000
    resources:
        mem_mb = 8000
    shell:
        "seq_depth=$(cat {input.seqdepth}); "
        "scale_factor=$(echo \"scale = 6; {params.scale_constant} / $seq_depth\" | bc); "
        "echo $scale_factor; "
        "bedtools genomecov -bg -i {input.bed} -scale $scale_factor -g {config[chrom_sizes]} > {output}"

rule bedgraph_to_bigwig:
    input: "{directory}{filename}.bedgraph"
    output: "{directory}{filename}.bigwig"
    threads: 1
    shell:
        "bedGraphToBigWig {input} {config[chrom_sizes]} {output}"

rule filter_excluded_regions:
    input: "{directory}{filename}.bed"
    output: "{directory}{filename}.filtered.bed"
    shell:
        "bedtools intersect -v -a {input} -b {config[exclusion_list]} > {output}"


rule filter_excluded_regions_narrowpeak:
    input: "{directory}{filename}.narrowPeak"
    output: "{directory}{filename}.filtered.narrowPeak"
    shell:
        "bedtools intersect -v -a {input} -b {config[exclusion_list]} > {output}"

rule call_seacr_peaks:
    input: OUTPUT_DIR + "{dir_type}/signal/{sample}.scaled.bedgraph"
    output: OUTPUT_DIR + "{dir_type}/peaks/{sample}.{stringency}.bed"
    threads: 1
    log: OUTPUT_DIR + "{dir_type}/peaks/{sample}.vs-ctrl.{stringency}.log"
    params:
        prefix = OUTPUT_DIR + "{dir_type}/peaks/{sample}.{stringency}",
        tmp = OUTPUT_DIR + "{dir_type}/peaks/{sample}.{stringency}.{stringency}.bed"
    shell:
        "SEACR_1.3.sh {input} 0.01 non {wildcards.stringency} {params.prefix} &> {log}; "
        "mv {params.tmp} {output}"


def get_expt_and_ctrl_bedgraphs(wildcards):
    if wildcards.dir_type == "samples":
        control_sample = CONTROLS[wildcards.sample]
    if wildcards.dir_type == "conditions":
        control_sample = wildcards.sample + "_CONTROL"
    expt_bg = os.path.join(OUTPUT_DIR, wildcards.dir_type, "signal", wildcards.sample + ".scaled.bedgraph")
    ctrl_bg = os.path.join(OUTPUT_DIR, wildcards.dir_type, "signal", control_sample + ".scaled.bedgraph")
    return {'expt' : expt_bg, 'ctrl' : ctrl_bg}

rule call_seacr_peaks_vs_control:
    input: unpack(get_expt_and_ctrl_bedgraphs)
    output: OUTPUT_DIR + "{dir_type}/peaks/{sample}.vs-ctrl.{stringency}.bed"
    threads: 1
    log: OUTPUT_DIR + "{dir_type}/peaks/{sample}.vs-ctrl.{stringency}.log"
    params:
        prefix = OUTPUT_DIR + "{dir_type}/peaks/{sample}.vs-ctrl.{stringency}",
        tmp = OUTPUT_DIR + "{dir_type}/peaks/{sample}.vs-ctrl.{stringency}.{stringency}.bed"
    shell:
        "SEACR_1.3.sh {input.expt} {input.ctrl} non {wildcards.stringency} {params.prefix} &> {log}; "
        "mv {params.tmp} {output}"


def get_expt_and_ctrl_bams(wildcards):
    if wildcards.dir_type == "samples":
        control_sample = CONTROLS[wildcards.sample]
    if wildcards.dir_type == "conditions":
        control_sample = wildcards.sample + "_CONTROL"
    expt_bg = os.path.join(OUTPUT_DIR, wildcards.dir_type, "align", wildcards.sample + ".cleaned.bam")
    ctrl_bg = os.path.join(OUTPUT_DIR, wildcards.dir_type, "align", control_sample + ".cleaned.bam")
    return {'expt' : expt_bg, 'ctrl' : ctrl_bg}

rule call_macs2_peaks_vs_control:
    input: unpack(get_expt_and_ctrl_bams)
    output: OUTPUT_DIR + "{dir_type}/peaks_macs2/{sample}.macs2_q0.1_peaks.narrowPeak"
    threads: 1
    log: OUTPUT_DIR + "{dir_type}/peaks_macs2/{sample}.macs2_q0.1.log"
    params:
        outdir = OUTPUT_DIR + "{dir_type}/peaks_macs2"
    shell:
        "macs2 callpeak -t {input.expt} -c {input.ctrl} "
        "-n {wildcards.sample}.macs2_q0.1 --outdir {params.outdir} "
        "-g hs -f BAMPE -q 0.1 &> {log}"


def get_samples_per_condition(wildcards):
    return [os.path.join(OUTPUT_DIR, "samples/align/", s + ".cleaned.bam") for s in CONDITIONS[wildcards.condition]]

rule merge_samples:
    input: get_samples_per_condition
    output: OUTPUT_DIR + "conditions/align/{condition}.cleaned.bam"
    threads: THREADS
    shell:
        # "cat {input} | sort -k1,1 -k2,2n -k3,3n > {output}"
        "samtools merge -@ {threads} {output} {input}"

def get_controls_per_condition(wildcards):
    return [os.path.join(OUTPUT_DIR, "samples/align/", s + ".cleaned.bam") for s in CONDITION_CONTROLS[wildcards.condition]]

rule merge_controls:
    input: get_controls_per_condition
    output: OUTPUT_DIR + "conditions/align/{condition}_CONTROL.cleaned.bam"
    threads: THREADS
    shell:
        # "cat {input} | sort -k1,1 -k2,2n -k3,3n > {output}"
        "samtools merge -@ {threads} {output} {input}"

def get_sample_spikeins_per_condition(wildcards):
    return [os.path.join(OUTPUT_DIR, "samples/align-spikein/", s + ".spikein.cleaned.bam") for s in CONDITIONS[wildcards.condition]]

rule merge_sample_spikeins:
    input: get_sample_spikeins_per_condition
    output: OUTPUT_DIR + "conditions/align-spikein/{condition}.spikein.cleaned.bam"
    threads: THREADS
    shell:
        # "cat {input} | sort -k1,1 -k2,2n -k3,3n > {output}"
        "samtools merge -@ {threads} {output} {input}"



def get_control_spikeins_per_condition(wildcards):
    return [os.path.join(OUTPUT_DIR, "samples/align-spikein/", s + ".spikein.cleaned.bam") for s in CONDITION_CONTROLS[wildcards.condition]]

rule merge_control_spikeins:
    input: get_control_spikeins_per_condition
    output: OUTPUT_DIR + "conditions/align-spikein/{condition}_CONTROL.spikein.cleaned.bam"
    threads: THREADS
    shell:
        # "cat {input} | sort -k1,1 -k2,2n -k3,3n > {output}"
        "samtools merge -@ {threads} {output} {input}"


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
