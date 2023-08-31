##=== R command ===## 
## Path to the project and histone list
# projPath = "/fh/fast/gottardo_r/yezheng_working/cuttag/CUTTag_tutorial"
# sampleList = c("K27me3_rep1", "K27me3_rep2", "K4me3_rep1", "K4me3_rep2", "IgG_rep1", "IgG_rep2")
# histList = c("K27me3", "K4me3", "IgG")

log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "message")
sink(log, type = "output")

library(tidyverse)

samples <- read_csv(snakemake@config$samplesheet)

alignment_summary <- map_dfr(samples$sample, function(s) {
	bowtie2_output <- file.path(snakemake@config$output_dir, "samples/align/", paste0(s, ".bowtie2.log"))
	align_results <- readLines(bowtie2_output)
	first_line <- which(str_detect(align_results, "reads; of these:"))[1]
	total_reads <- as.integer(str_extract(align_results[first_line], "^\\s*\\d+"))
	unaligned_reads <- as.integer(str_extract(align_results[first_line + 2], "^\\s*\\d+"))
	aligned_once <- as.integer(str_extract(align_results[first_line + 3], "^\\s*\\d+"))
	aligned_multi <- as.integer(str_extract(align_results[first_line + 4], "^\\s*\\d+"))
	total_aligned <- aligned_once + aligned_multi

	return(tibble(sample = s, total_reads = total_reads, aligned_total = total_aligned))

})

alignment_summary <- alignment_summary %>%
		mutate(aligned_unique = map_int(sample,
			~ as.integer(readLines(file.path(snakemake@config$output_dir, "samples/align", paste0(., ".cleaned.bam.seqdepth")))[1])))

if(!is.null(snakemake@config$reference_spikein) && snakemake@config$reference_spikein != "") {
	alignment_summary <- alignment_summary %>%
		mutate(spikein_aligned = map_int(sample, function(s) {
			bowtie2_output <- file.path(snakemake@config$output_dir, "samples/align-spikein/", paste0(s, ".spikein.bowtie2.log"))
			align_results <- readLines(bowtie2_output)
			first_line <- which(str_detect(align_results, "reads; of these:"))[1]
			aligned_once <- as.integer(str_extract(align_results[first_line + 3], "^\\s*\\d+"))
			aligned_multi <- as.integer(str_extract(align_results[first_line + 4], "^\\s*\\d+"))
			total_aligned <- aligned_once + aligned_multi
			return(total_aligned)
			})) %>%
		mutate(spikein_unique = map_int(sample,
			~ as.integer(readLines(file.path(snakemake@config$output_dir, "samples/align-spikein", paste0(., ".spikein.cleaned.bam.seqdepth")))[1])))

}

write_csv(alignment_summary, snakemake@output$alignment_summary)



# ## Collect the alignment results from the bowtie2 alignment summary files
# alignResult = c()
# for(hist in sampleList){
#   alignRes = read.table(paste0(projPath, "/alignment/sam/bowtie2_summary/", hist, "_bowtie2.txt"), header = FALSE, fill = TRUE)
#   alignRate = substr(alignRes$V1[6], 1, nchar(as.character(alignRes$V1[6]))-1)
#   histInfo = strsplit(hist, "_")[[1]]
#   alignResult = data.frame(Histone = histInfo[1], Replicate = histInfo[2], 
#                            SequencingDepth = alignRes$V1[1] %>% as.character %>% as.numeric, 
#                            MappedFragNum_hg38 = alignRes$V1[4] %>% as.character %>% as.numeric + alignRes$V1[5] %>% as.character %>% as.numeric, 
#                            AlignmentRate_hg38 = alignRate %>% as.numeric)  %>% rbind(alignResult, .)
# }
# alignResult$Histone = factor(alignResult$Histone, levels = histList)
# alignResult %>% mutate(AlignmentRate_hg38 = paste0(AlignmentRate_hg38, "%"))

# spikeAlign = c()
# for(hist in sampleList){
#   spikeRes = read.table(paste0(projPath, "/alignment/sam/bowtie2_summary/", hist, "_bowtie2_spikeIn.txt"), header = FALSE, fill = TRUE)
#   alignRate = substr(spikeRes$V1[6], 1, nchar(as.character(spikeRes$V1[6]))-1)
#   histInfo = strsplit(hist, "_")[[1]]
#   spikeAlign = data.frame(Histone = histInfo[1], Replicate = histInfo[2], 
#                           SequencingDepth = spikeRes$V1[1] %>% as.character %>% as.numeric, 
#                           MappedFragNum_spikeIn = spikeRes$V1[4] %>% as.character %>% as.numeric + spikeRes$V1[5] %>% as.character %>% as.numeric, 
#                           AlignmentRate_spikeIn = alignRate %>% as.numeric)  %>% rbind(spikeAlign, .)
# }
# spikeAlign$Histone = factor(spikeAlign$Histone, levels = histList)
# spikeAlign %>% mutate(AlignmentRate_spikeIn = paste0(AlignmentRate_spikeIn, "%"))