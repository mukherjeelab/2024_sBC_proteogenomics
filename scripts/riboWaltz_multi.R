#!/usr/bin/env Rscript
options(echo = TRUE)

library(riboWaltz)
library(here)
library(tidyverse)

sink("output.log", append = TRUE)

#Make the annotation file once, then can just load instead. 
#annotation_file <- "/beevol/home/bartholk/Ribotricer/gtf/gencode.v26.primary_assembly.annotation.gtf"
#annot <- create_annotation(annotation_file)
#saveRDS(annot, file = "annot_object.rds") # Save the annotation object



#load annotation already made.
annot <- readRDS(here("annot_object.rds"))

#created this large rds object and saved. Took ~20 min. Don't want to rerun. 
#reads_list <- bamtolist(bamfolder = "bam/", annotation=annot)
#saveRDS(reads_list, file = "2024_RiboALLSamples.rds")


#setup files needed for most downstream functions
reads_list <- readRDS(here("2024_RiboALLSamples.rds"))
psite_offset <- psite(reads_list, flanking = 6, extremity = "auto")
reads_psite_list <- psite_info(reads_list, psite_offset)

rlistDups <- lapply(reads_psite_list, function(x) x[!duplicated(x), ])
print(rlistDups)

#codon_coverage <- codon_coverage(reads_psite_list, annot, psite = TRUE)
#write_tsv(codon_coverage, "summary/codon_coverage_ALL.tsv.gz")


#cds_coverage <- cds_coverage(reads_psite_list, annot)
#write_tsv(cds_coverage, "summary/cds_coverage_ALL.tsv.gz")


#first look at frag length of all samples independently
input_samples_1 <- c("NM2022_0001_RP_iPSC", "NM2022_0002_RP_iPSC", "NM2022_0003_RP_DE", "NM2022_0004_RP_DE")
input_samples_2 <- c("NM2022_0005_IP_iPSC", "NM2022_0006_IP_iPSC", "NM2022_0007_IP_DE", "NM2022_0008_IP_DE", "NM2022_0085_IP_beta", "NM2022_0086_IP_beta")

example_length_dist_1 <- rlength_distr(reads_list, sample = input_samples_1, multisamples = "independent", plot_style = "facet", cl = 99)
length_dist_plt_1 <- example_length_dist_1[["plot"]] + theme_classic()+ theme(legend.position = "none")
ggsave(plot = length_dist_plt_1, file = "plots/average_length_dist_ALL_1.pdf")

example_length_dist_2 <- rlength_distr(reads_list, sample = input_samples_2, multisamples = "independent", plot_style = "facet", cl = 99)
length_dist_plt_2 <- example_length_dist_2[["plot"]] + theme_classic() + theme(legend.position = "none")
ggsave(plot = length_dist_plt_2, file = "plots/average_length_dist_ALL_2.pdf")






#now I want to group. So first group is just comparing methods (RP vs IP) for each timepoint. I am averaging replicates.
input_samples_iPSC <- list("RP_iPSC" = c("NM2022_0001_RP_iPSC", "NM2022_0002_RP_iPSC"),
  			"IP_iPSC" = c("NM2022_0005_IP_iPSC", "NM2022_0006_IP_iPSC"))
example_length_dist_iPSC <- rlength_distr(reads_list, sample = input_samples_iPSC, multisamples = "average", plot_style = "dodge", cl = 99)
length_dist_plt_iPSC <- example_length_dist_iPSC[["plot"]] + theme_classic()
ggsave(plot = length_dist_plt_iPSC, file = "plots/average_length_dist_iPSC.pdf")


input_samples_DE <- list("RP_DE" = c("NM2022_0003_RP_DE", "NM2022_0004_RP_DE"),
  			"IP_DE" = c("NM2022_0007_IP_DE", "NM2022_0008_IP_DE"))
example_length_dist_DE <- rlength_distr(reads_list, sample = input_samples_DE, multisamples = "average", plot_style = "dodge", cl = 99)
length_dist_plt_DE <- example_length_dist_DE[["plot"]] + theme_classic()
ggsave(plot = length_dist_plt_DE, file = "plots/average_length_dist_DE.pdf")


#now I am looking at all the samples for IP method, across all cell types. Again, replicates averaged. 
input_samples_IPmethod <- list("IP_iPSC" = c("NM2022_0005_IP_iPSC", "NM2022_0006_IP_iPSC"),
  			"IP_DE" = c("NM2022_0007_IP_DE", "NM2022_0008_IP_DE"), 
  			"IP_Beta" = c("NM2022_0085_IP_beta", "NM2022_0086_IP_beta"))

example_length_dist_IP <- rlength_distr(reads_list, sample = input_samples_IPmethod, multisamples = "average", plot_style = "dodge", cl = 99)
length_dist_plt_IP <- example_length_dist_IP[["plot"]] + theme_classic() 
ggsave(plot = length_dist_plt_IP, file = "plots/average_length_dist_IPmethod.pdf")




#Next plot type. Looking at the regions that the psites map to. Doing the same thing where I do each indepenently and then look at some replicate averages/groups
input_samples <- c("NM2022_0001_RP_iPSC", "NM2022_0002_RP_iPSC","NM2022_0003_RP_DE", "NM2022_0004_RP_DE","NM2022_0005_IP_iPSC", "NM2022_0006_IP_iPSC", "NM2022_0007_IP_DE", "NM2022_0008_IP_DE", "NM2022_0085_IP_beta", "NM2022_0086_IP_beta")

psite_per_region <- region_psite(reads_psite_list, annot, sample = input_samples, multisamples = "independent", plot_style = "stack", cl = 85)
psite_per_region_plt <- psite_per_region[["plot"]] + theme_classic()
ggsave(plot = psite_per_region_plt, file = "plots/psite_per_region_ALL.pdf")


psite_per_region_iPSC <- region_psite(reads_psite_list, annot, sample = input_samples_iPSC, multisamples = "average", plot_style = "dodge", cl = 85)
psite_per_region_plt_iPSC <- psite_per_region_iPSC[["plot"]] + theme_classic()
ggsave(plot = psite_per_region_plt_iPSC, file = "plots/psite_per_region_iPSC.pdf")

psite_per_region_DE <- region_psite(reads_psite_list, annot, sample = input_samples_DE, multisamples = "average", plot_style = "dodge", cl = 85)
psite_per_region_plt_DE <- psite_per_region_DE[["plot"]] + theme_classic()
ggsave(plot = psite_per_region_plt_DE, file = "plots/psite_per_region_DE.pdf")

psite_per_region_IP <- region_psite(reads_psite_list, annot, sample = input_samples_IPmethod, multisamples = "average", plot_style = "dodge", cl = 85)
psite_per_region_plt_IP <- psite_per_region_IP[["plot"]] + theme_classic()
ggsave(plot = psite_per_region_plt_IP, file = "plots/psite_per_region_IPmethod.pdf")



frames_psites <- frame_psite(reads_psite_list, annot, sample = input_samples, multisamples = "independent", plot_style = "facet", region = "cds")
frames_psites_plt <- frames_psites[["plot"]] + theme_classic()
ggsave(plot = frames_psites_plt, file = "plots/frames_psites_ALL.pdf")

frames_psites_iPSC <- frame_psite(reads_psite_list, annot, sample = input_samples_iPSC, multisamples = "average", plot_style = "dodge", region = "cds")
frames_psites_iPSC_plt <- frames_psites_iPSC[["plot"]] + theme_classic()
ggsave(plot = frames_psites_iPSC_plt, file = "plots/frames_psites_iPSC.pdf")

frames_psites_DE <- frame_psite(reads_psite_list, annot, sample = input_samples_DE, multisamples = "average", plot_style = "dodge", region = "cds")
frames_psites_DE_plt <- frames_psites_DE[["plot"]] + theme_classic()
ggsave(plot = frames_psites_DE_plt, file = "plots/frames_psites_DE.pdf")

frames_psites_IP <- frame_psite(reads_psite_list, annot, sample = input_samples_IPmethod, "average", plot_style = "dodge", region = "cds")
frames_psites_IP_plt <- frames_psites_IP[["plot"]] + theme_classic()
ggsave(plot = frames_psites_IP_plt, file = "plots/frames_psites_IPmethod.pdf")




metaprofile <- metaprofile_psite(reads_psite_list, annot, sample = input_samples, multisamples = "independent", plot_style = "facet", utr5l = 20, cdsl = 40, utr3l = 20)
metaprofile_plt <- metaprofile[["plot"]] + theme_classic()
ggsave(plot = metaprofile_plt, file = "plots/metaprofile_ALL.pdf")


metaprofile_iPSC <- metaprofile_psite(reads_psite_list, annot, sample = input_samples_iPSC, multisamples = "average", plot_style = "overlap", utr5l = 20, cdsl = 40, utr3l = 20)
metaprofile_iPSC_plt <- metaprofile_iPSC[["plot"]] + theme_classic()
ggsave(plot = metaprofile_iPSC_plt, file = "plots/metaprofile_iPSC.pdf")

metaprofile_DE <- metaprofile_psite(reads_psite_list, annot, sample = input_samples_DE, multisamples = "average", plot_style = "overlap", utr5l = 20, cdsl = 40, utr3l = 20)
metaprofile_DE_plt <- metaprofile_DE[["plot"]] + theme_classic()
ggsave(plot = metaprofile_DE_plt, file = "plots/metaprofile_DE.pdf")

metaprofile_IP <- metaprofile_psite(reads_psite_list, annot, sample = input_samples_IPmethod, multisamples = "average", plot_style = "overlap", utr5l = 20, cdsl = 40, utr3l = 20)
metaprofile_IP_plt <- metaprofile_IP[["plot"]] + theme_classic()
ggsave(plot = metaprofile_IP_plt, file = "plots/metaprofile_IPmethod.pdf")



metaheatmap <- metaheatmap_psite(reads_psite_list, annot, sample = input_samples, multisamples = "independent", utr5l = 20, cdsl = 40, utr3l = 20)
metaheatmap_plt <- metaheatmap[["plot"]] + theme_classic()
ggsave(plot = metaheatmap_plt, file = "plots/metaheatmap_ALL.pdf")

metaheatmap_iPSC <- metaheatmap_psite(reads_psite_list, annot, sample = input_samples_iPSC, multisamples = "average", utr5l = 20, cdsl = 40, utr3l = 20)
metaheatmap_iPSC_plt <- metaheatmap_iPSC[["plot"]] + theme_classic()
ggsave(plot = metaheatmap_iPSC_plt, file = "plots/metaheatmap_iPSC.pdf")

metaheatmap_DE <- metaheatmap_psite(reads_psite_list, annot, sample = input_samples_DE, multisamples = "average", utr5l = 20, cdsl = 40, utr3l = 20)
metaheatmap_DE_plt <- metaheatmap_DE[["plot"]] + theme_classic()
ggsave(plot = metaheatmap_DE_plt, file = "plots/metaheatmap_DE.pdf")

metaheatmap_IP <- metaheatmap_psite(reads_psite_list, annot, sample = input_samples_IPmethod, multisamples = "average", utr5l = 20, cdsl = 40, utr3l = 20)
metaheatmap_IP_plt <- metaheatmap_IP[["plot"]] + theme_classic()
ggsave(plot = metaheatmap_IP_plt, file = "plots/metaheatmap_IPmethod.pdf")
sink()
