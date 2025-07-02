#Load packages and functions----
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(BiocParallel)
library(progress)
library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
args <- commandArgs(trailingOnly = TRUE)
script_folder <- args[1]
params_file <- args[2]
source(file.path(script_folder, "CrossPeak_functions.R"))
source(file.path(params_file))
print(paste0('Using parameters file ', params_file))

#Load necessary rds files from previous steps----
species_spec_peaks <- readRDS(paste0(folder,'species_spec_peaks.rds'))
species_spec_summits <- readRDS(paste0(folder,'species_spec_summits.rds'))

#Step 7: Categorizing putative species-specific peaks into peaks that are truly species-specific vs those that can be added to consensus list----
indel_values <- c(round2_indel_high = indel_high, round2_indel_med = indel_med,round2_indel_low = indel_low, round2_summit_deletion = indel_high)
result = getRound2PeakCategories(species,species_spec_peaks, species_spec_summits, indel_low, indel_med, indel_high,
                                 genome_names,folder, peak_half_width, summit_half_width_rd2) #see function description for details of categories. category metadata will be added to final peak sets
species_consensus_categories1 = result$species_consensus_categories #List of lists: $species1$species2 will be species2 coordinates for species1 peaks that were lifted to species2 (in a consensus category)
species_specific_categories = result$species_specific_categories #List of lists: $species1$species2 will be species1 coordinates for species1 peaks that failed to lift to species2 (or lifted but with large indels such that the peak cannot be confidently located)
#Export putative round 2 consensus peaks
for (sp in species){
  other_species = setdiff(species, sp)
  for (sp_other in other_species){
    export.bed(species_consensus_categories1[[sp]][[sp_other]],paste0(folder,sp,'_on_',sp_other,'_rd2_con_peaks.bed'))
  }
}
saveRDS(species_consensus_categories1, paste0(folder,'species_consensus_categories1.rds'))
saveRDS(species_specific_categories, paste0(folder,'species_specific_categories.rds'))