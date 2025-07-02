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
summits <- readRDS(paste0(folder,'summits.rds'))
peaks <- readRDS(paste0(folder,'peaks.rds'))
species_blacklists <- readRDS(paste0(folder,'species_blacklists.rds'))
consensus_summits <- readRDS(paste0(folder,'consensus_summits.rds'))
consensus_peaks1 <-  readRDS(paste0(folder,'consensus_peaks1.rds'))
species_merged_metadata <- readRDS(paste0(folder,'species_merged_metadata.rds'))
species_on_ancestral_peaks <- readRDS(paste0(folder,'species_on_ancestral_peaks.rds'))

#Step 4: Import halliftover results of consensus peaks back to original species genomes, concatenate summit regions, and make peaks----
#Input filenames should be [species]_consensus_summits_lo.bed
species_con_summits1 <- vector(mode = 'list', length = length(species)); names(species_con_summits1) <- species
species_con_summits2 <- vector(mode = 'list', length = length(species)); names(species_con_summits2) <- species
species_con_summits <- vector(mode = 'list', length = length(species)); names(species_con_summits) <- species
species_con_summits_failed <- vector(mode = 'list', length = length(species)); names(species_con_summits_failed) <- species
species_rd1_peaks1 <- vector(mode = 'list', length = length(species)); names(species_rd1_peaks1) <- species
problem_summits  <- vector(mode = 'list', length = length(species)); names(problem_summits) <- species
for (sp in species){
  print(paste0('Processing ', sp, ' peaks:'))
  result = concatenateAndSummarizeLiftovers(consensus_summits,paste0(folder,sp,'_consensus_summits'), summit_half_width, indel_tol)
  species_con_summits1[[sp]] = result[[1]]; genome(species_con_summits1[[sp]]) = genome_names[[sp]]
  species_con_summits2[[sp]] = result[[2]]; genome(species_con_summits2[[sp]]) = genome_names[[sp]]
  species_con_summits_failed[[sp]] = result[[3]]
  species_con_summits[[sp]] <- sortGrByName(species_con_summits2[[sp]])
  #Add info about peak origin to species_rd1_peaks (is not retained in halliftover on command line)
  matching_indices <- match(names(species_con_summits[[sp]]), names(consensus_summits))
  names_to_add <- consensus_summits$orig_name[matching_indices]
  species_con_summits[[sp]]$orig_name <- names_to_add
  #Make peaks from summits
  species_rd1_peaks1[[sp]] <- makePeaksFromSummitRanges(species_con_summits[[sp]],peak_half_width)
  #Find peaks that do not map back to the original species peak (only considering region around summit),
  #or are too far from the original peak (> summit_adjust/10*indel_tol)
  problem_summits[[sp]] = compareMultiConsensusAndOriginalSummits(summits[[sp]],species_merged_metadata[[sp]],species_con_summits[[sp]], indel_tol, summit_half_width, prefixes[[sp]])
} 
saveRDS(species_con_summits1,paste0(folder,'species_con_summits1.rds'))
saveRDS(species_con_summits2,paste0(folder,'species_con_summits2.rds'))
saveRDS(species_con_summits,paste0(folder,'species_con_summits.rds'))
saveRDS(species_con_summits_failed,paste0(folder,'species_con_summits_failed.rds'))
saveRDS(species_rd1_peaks1, paste0(folder,'species_rd1_peaks1.rds'))

#Step 5: Exclude peaks that were too far from the original location or overlap blacklisted regions in either species----
result = excludeFailedProblemBlacklistPeaks(consensus_peaks1,species_rd1_peaks1,problem_summits, species_blacklists,species, remove_non_standard_chr, chr_to_remove)
consensus_peaks2 = result[[1]]; species_rd1_peaks2 = result[[2]]
consensus_peaks = consensus_peaks2;
consensus_peaks_success = consensus_peaks[consensus_peaks$liftback_result=='',]
saveRDS(consensus_peaks,paste0(folder,'consensus_peaks.rds'))
saveRDS(species_rd1_peaks2,paste0(folder,'species_rd1_peaks2.rds'))

#Sort by name so all species are in same order
species_rd1_peaks = species_rd1_peaks2
for (sp in species){
  species_rd1_peaks[[sp]] = sortGrByName(species_rd1_peaks[[sp]])
  strand(species_rd1_peaks[[sp]]) <-  '*' #not considering strand info for quantification
}
saveRDS(species_rd1_peaks,paste0(output_folder,'species_rd1_peaks.rds'))

#Export bed files if user wishes to end the pipeline after step 5
for (sp in species){
  export.bed(species_rd1_peaks[[sp]],paste0(output_folder,sp,'_rd1_peaks.bed'))
}

#Make consensus data frame for faster computation in next steps
consensus_df <- as.data.frame(consensus_peaks)
consensus_df <- consensus_df %>%
  separate_rows(orig_name, orig_species, sep = ";") %>%
  separate_rows(liftback_result, liftback_species, sep = ";")
saveRDS(consensus_df,paste0(folder,'consensus_df.rds'))

#Step 6: combining species-specific peaks that failed to make it into consensus set----
#Use larger summit half width to increase chances of successful liftover and variable tolerances for different categories
species_spec_peaks <- getMultiSpeciesSpecificPeaks(species, species_on_ancestral_peaks, peaks, consensus_peaks, consensus_df)
saveRDS(species_spec_peaks, paste0(folder,'species_spec_peaks.rds'))

#Export putative species-specific summits and peaks
species_spec_summits  <- vector(mode = 'list', length = length(species)); names(species_spec_summits) <- species
for (sp in species){
  export.bed(species_spec_peaks[[sp]],paste0(folder,sp,'_spec_peaks.bed'))
  species_spec_summits[[sp]] <- makeSummitsFromPeaks(species_spec_peaks[[sp]], peak_half_width,summit_half_width_rd2)
  export.bed(species_spec_summits[[sp]],paste0(folder,sp,'_spec_summits.bed'))
}
saveRDS(species_spec_summits, paste0(folder,'species_spec_summits.rds'))