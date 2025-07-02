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

#Step 2: Import halliftover results on ancestral genome, concatenate summit regions, and make peaks----
#Input filenames should be same as original, but ending in "lo"
species_on_ancestral_summits1 <- vector(mode = 'list', length = length(species)); names(species_on_ancestral_summits1) <- species
species_on_ancestral_summits2 <- vector(mode = 'list', length = length(species)); names(species_on_ancestral_summits2) <- species
species_on_ancestral_summits <- vector(mode = 'list', length = length(species)); names(species_on_ancestral_summits) <- species
species_on_ancestral_summits_failed <- vector(mode = 'list', length = length(species)); names(species_on_ancestral_summits_failed) <- species
species_on_ancestral_peaks1 <- vector(mode = 'list', length = length(species)); names(species_on_ancestral_peaks1) <- species
species_on_ancestral_peaks2 <- vector(mode = 'list', length = length(species)); names(species_on_ancestral_peaks2) <- species
species_on_ancestral_peaks <- vector(mode = 'list', length = length(species)); names(species_on_ancestral_peaks) <- species
for (sp in species){
  #Import and summarize
  print(paste0('Processing ', sp, ' peaks:'))
  result <- concatenateAndSummarizeLiftovers(summits[[sp]],paste0(folder,sp,'_summits'), summit_half_width, indel_tol)
  species_on_ancestral_summits1[[sp]] = result[[1]]; genome(species_on_ancestral_summits1[[sp]]) = genome_ancestral
  species_on_ancestral_summits2[[sp]] = result[[2]]; genome(species_on_ancestral_summits2[[sp]]) = genome_ancestral
  species_on_ancestral_summits_failed[[sp]] = result[[3]]
  #Make peaks from summits
  species_on_ancestral_peaks1[[sp]] <- makePeaksFromSummitRanges(species_on_ancestral_summits2[[sp]],peak_half_width)
  #Exclude self-overlapping peaks (ie peaks that mapped to same region in ancestor)
  species_on_ancestral_peaks2[[sp]] <- excludeSelfOverlaps(species_on_ancestral_peaks1[[sp]],2*indel_tol, ignore_strand = TRUE)
  species_on_ancestral_peaks[[sp]] <- species_on_ancestral_peaks2[[sp]]
}
saveRDS(species_on_ancestral_summits1, paste0(folder,'species_on_ancestral_summits1.rds'))
saveRDS(species_on_ancestral_summits2, paste0(folder,'species_on_ancestral_summits2.rds'))
saveRDS(species_on_ancestral_summits_failed, paste0(folder,'species_on_ancestral_summits_failed.rds'))
saveRDS(species_on_ancestral_summits, paste0(folder,'species_on_ancestral_summits.rds'))
saveRDS(species_on_ancestral_peaks1, paste0(folder,'species_on_ancestral_peaks1.rds'))
saveRDS(species_on_ancestral_peaks2, paste0(folder,'species_on_ancestral_peaks2.rds'))
saveRDS(species_on_ancestral_peaks, paste0(folder,'species_on_ancestral_peaks.rds'))

#Check overlaps (optional - can help in deciding overlap_tol)
plot_overlap_histograms(species_on_ancestral_peaks, species, save_plots = TRUE, figure_folder)

# Step 3: Merging peaks on the ancestral genome----
result <- mergePeaksAcrossMultipleSpecies(species_on_ancestral_peaks, species, prefixes, overlap_tol, peak_half_width)
consensus_peaks1 = result$mergedPeaks; species_merged_metadata = result$metadata
genome(consensus_peaks1) = genome_ancestral
saveRDS(consensus_peaks1,paste0(folder,'consensus_peaks1','.rds'))
saveRDS(species_merged_metadata, paste0(folder,'species_merged_metadata', '.rds'))
#Make summits and export for liftover
consensus_summits <- makeSummitsFromPeaks(consensus_peaks1, peak_half_width, summit_half_width)
saveRDS(consensus_summits,paste0(folder,'consensus_summits','.rds'))
export.bed(consensus_summits, paste0(folder,'consensus_summits','.bed')) #note that due to indexing differences, start position will be -1 in bed file (but will be corrected when imported back to R)

