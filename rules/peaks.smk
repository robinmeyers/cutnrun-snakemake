
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
    resources:
        mem_mb = 8000
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
    resources:
        mem_mb = 8000
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
    resources:
        mem_mb = 8000
    shell:
        "macs2 callpeak -t {input.expt} -c {input.ctrl} "
        "-n {wildcards.sample}.macs2_q0.1 --outdir {params.outdir} "
        "-g hs -f BAMPE -q 0.1 &> {log}"

