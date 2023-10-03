
rule align:
    input:
        r1 = OUTPUT_DIR + "samples/trim/{sample}_R1.fastq.gz",
        r2 = OUTPUT_DIR + "samples/trim/{sample}_R2.fastq.gz"
    output:
        sam = temp(OUTPUT_DIR + "samples/align/{sample}.sam")
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
        sam = temp(OUTPUT_DIR + "samples/align-spikein/{sample}.spikein.sam")
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




rule filter_excluded_regions_bam:
    input:
        OUTPUT_DIR + "samples/align/{sample}.bam"
    output:
        temp(OUTPUT_DIR + "samples/align/{sample}.filtered.bam")
    params:
    threads: THREADS
    resources:
        # mem_mb = 32000
    shell:
        "bedtools intersect -v -a {input} -b {config[exclusion_list]} | "
        "samtools fixmate -@ {threads} -m - - | "
        "samtools view -@ {threads} -b -f 0x3 -F 0x4 - > {output}"


rule sort_filter_bam:
    input:
        OUTPUT_DIR + "samples/align/{sample}.filtered.bam"
    output:
        OUTPUT_DIR + "samples/align/{sample}.cleaned.bam"
    params:
        tmp = OUTPUT_DIR + "samples/align/{sample}"
    threads: THREADS
    resources:
        mem_mb = 32000
    shell:
        # "samtools fixmate -m -@ {threads} {input} - | "
        "samtools sort -@ {threads} -T {params.tmp} {input} | "
        "samtools markdup - - | "
        "samtools view -b -f 0x3 -F 0x400 - > {output}"
        # "samtools sort -n - > {output}"
        


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
        "samtools view -b -f 0x3 -F 0x400 - > {output}"
        # "samtools sort -n - > {output}"

rule bam_to_bed:
    input:
        bam = OUTPUT_DIR + "{dir_type}/align/{sample}.cleaned.bam",
        bai = OUTPUT_DIR + "{dir_type}/align/{sample}.cleaned.bam.bai"
    output: temp(OUTPUT_DIR + "{dir_type}/align/{sample}.bed")
    threads: 1
    shell:
        "samtools collate -f -O {input.bam} | bedtools bamtobed -bedpe > {output}"

rule clean_bed:
    input: OUTPUT_DIR + "{dir_type}/align/{sample}.bed"
    output: OUTPUT_DIR + "{dir_type}/align/{sample}.cleaned.bed.gz"
    threads: 1
    shell:
        "awk '$1==$4 && $6-$2 < 1000 {{print $0}}' {input} | "
        "cut -f 1,2,6,7,8,9 | sort -k1,1 -k2,2n -k3,3n | "
        "gzip -c > {output} "



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