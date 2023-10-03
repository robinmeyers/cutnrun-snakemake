
rule bed_to_bedgraph:
    input: OUTPUT_DIR + "{dir_type}/align/{sample}.cleaned.bed.gz"
    output: temp(OUTPUT_DIR + "{dir_type}/signal/{sample}.bedgraph")
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
        bed = OUTPUT_DIR + "{dir_type}/align/{sample}.cleaned.bed.gz",
        seqdepth = OUTPUT_DIR + "{dir_type}/align/{sample}.cleaned.bam.seqdepth"
    output: OUTPUT_DIR + "{dir_type}/signal/{sample}.scaled.bedgraph"
    threads: 1
    params: scale_constant = 10000000
    resources:
        mem_mb = 8000
    shell:
        "seq_depth=$(cat {input.seqdepth}); "
        "scale_factor=$(echo \"scale = 6; {params.scale_constant} / $seq_depth\" | bc); "
        "echo $scale_factor; "
        "bedtools genomecov -bg -i {input.bed} -scale $scale_factor -g {config[chrom_sizes]} > {output}"


rule bed_to_scaled_bedgraph_spikein:
    input:
        bed = OUTPUT_DIR + "{dir_type}/align/{sample}.cleaned.bed.gz",
        seqdepth = OUTPUT_DIR + "{dir_type}/align-spikein/{sample}.spikein.cleaned.bam.seqdepth"
    output: OUTPUT_DIR + "{dir_type}/signal-spikein/{sample}.scaled.bedgraph"
    threads: 1
    params: scale_constant = 10000
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
