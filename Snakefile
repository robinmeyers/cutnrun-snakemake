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


include: "rules/reads.smk"
include: "rules/align.smk"
include: "rules/signal.smk"
include: "rules/peaks.smk"
include: "rules/conditions.smk"

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

    if config['reference_spikein']:
        targets = targets + expand("samples/signal-spikein/{sample}.scaled.bigwig", sample=SAMPLES.keys())

    targets = targets + expand("samples/align/{sample}.cleaned.fragmentsize.txt", sample = SAMPLES.keys())
    targets = targets + ["samples/alignment_summary.csv"]

    targets = [os.path.join(OUTPUT_DIR, t) for t in targets]
    return targets

rule all:
    input: get_target_files
    run:
        print("workflow complete!")



def get_alignment_summary_inputs(wildcards):
    targets = []
    targets = targets + expand(OUTPUT_DIR + "samples/align/{sample}.filtered.bam", sample=SAMPLES.keys())
    targets = targets + expand(OUTPUT_DIR + "samples/align/{sample}.cleaned.bam", sample=SAMPLES.keys())

    if "control" in samples.columns:
        targets = targets + expand(OUTPUT_DIR + "samples/peaks_macs2/{sample}.macs2_q0.1_peaks.narrowPeak", sample=CONTROLS.keys())

    if config['reference_spikein']:
        targets = targets + expand(OUTPUT_DIR + "samples/align-spikein/{sample}.spikein.cleaned.bam.seqdepth", sample=SAMPLES.keys())

    return targets

rule alignment_summary:
    input: get_alignment_summary_inputs
        
    output:
        alignment_summary = OUTPUT_DIR + "samples/alignment_summary.csv"
    log: OUTPUT_DIR + "samples/alignment_summary.log"
    script: "scripts/alignment_summary.R"

# rule mapping_stats:
#     input:
#     output:
#     shell:




