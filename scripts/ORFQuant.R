#!/usr/bin/env Rscript

library("ORFquant")
library("here")

#prepare_annotation_files(annotation_directory = ".", twobit_file = here("GRCh38.primary_assembly.genome.2bit"), gtf_file = here("gencode.v26.primary_assembly.annotation.gtf"), export_bed_tables_TxDb = F,forge_BSgenome = F,create_TxDb = T, genome_seq = "/beevol/home/neellab/genomes/human/hg38/fasta/GRCh38.primary_assembly.genome.fa")

run_ORFquant(for_ORFquant_file = here("Betacell_mergedBamsAligned.sortedByCoord.out.bam_for_SaTAnn"), annotation_file = "gencode.v26.primary_assembly.annotation.gtf_Rannot",n_cores = 6)
