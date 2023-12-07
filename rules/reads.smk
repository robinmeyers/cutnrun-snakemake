
def get_paired_fqs(wildcards):
    if (wildcards.sample not in SAMPLES.keys()):
        return {'r1' : [], 'r2': []}
        # raise ValueError(wildcards.sample + " is not a sample in the samplesheet")
    fastq_ids = SAMPLES[wildcards.sample]
    fastq_regex = "(" + "|".join(fastq_ids) + ")"

    r1_fastq_regex = ".*/" + fastq_regex + "(_S[0-9]+)?(_L[0-9]+)?_R1(_001)?.fastq.gz$"
    r2_fastq_regex = ".*/" + fastq_regex + "(_S[0-9]+)?(_L[0-9]+)?_R2(_001)?.fastq.gz$"
    
    fastq_dir_files = glob.glob(os.path.join(FASTQ_DIR, "*"), recursive=True)

    r1 = list(filter(re.compile(r1_fastq_regex).search, fastq_dir_files))
    r2 = list(filter(re.compile(r2_fastq_regex).search, fastq_dir_files))

    print(r1)
    print(r2)

    # r1 = glob.glob(os.path.join(FASTQ_DIR, "**", sample_id + "_*R1*.fastq.gz"),
    #     recursive=True)
    # r2 = glob.glob(os.path.join(FASTQ_DIR, "**", sample_id + "_*R2*.fastq.gz"), 
    #     recursive=True)
    if len(r1) == 0:
        raise ValueError(wildcards.sample + " has no matching input fastq file: " + fastq_regex)
    if len(r1) != len(r2):
        raise ValueError(wildcards.sample + " has different numbers of R1 and R2 fastq files: " + fastq_regex)
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
