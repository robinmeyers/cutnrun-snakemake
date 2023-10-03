
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
