
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
