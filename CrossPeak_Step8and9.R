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
species_blacklists <- readRDS(paste0(folder,'species_blacklists.rds'))
peaks_orig <- readRDS(paste0(folder,'peaks_orig.rds'))
peaks <- readRDS(paste0(folder,'peaks.rds'))
species_on_ancestral_summits1 <- readRDS(paste0(folder,'species_on_ancestral_summits1.rds'))
species_on_ancestral_summits2 <- readRDS(paste0(folder,'species_on_ancestral_summits2.rds'))
species_on_ancestral_peaks <- readRDS(paste0(folder,'species_on_ancestral_peaks.rds'))
species_on_ancestral_peaks1 <- readRDS(paste0(folder,'species_on_ancestral_peaks1.rds'))
species_on_ancestral_peaks2 <- readRDS(paste0(folder,'species_on_ancestral_peaks2.rds'))
consensus_peaks <- readRDS(paste0(folder,'consensus_peaks.rds'))
consensus_df <- readRDS(paste0(folder,'consensus_df.rds'))
species_spec_peaks <- readRDS(paste0(folder,'species_spec_peaks.rds'))
species_spec_summits <- readRDS(paste0(folder,'species_spec_summits.rds'))
species_con_summits1 <- readRDS(paste0(folder,'species_con_summits1.rds'))
species_con_summits2 <- readRDS(paste0(folder,'species_con_summits2.rds'))
species_rd1_peaks <- readRDS(paste0(output_folder,'species_rd1_peaks.rds'))
species_rd1_peaks1 <- readRDS(paste0(folder,'species_rd1_peaks1.rds'))
species_rd1_peaks2 <- readRDS(paste0(folder,'species_rd1_peaks2.rds'))
species_consensus_categories1 <- readRDS(paste0(folder,'species_consensus_categories1.rds'))
species_specific_categories <- readRDS( paste0(folder,'species_specific_categories.rds'))
indel_values <- c(round2_indel_high = indel_high, round2_indel_med = indel_med,round2_indel_low = indel_low, round2_summit_deletion = indel_high)

#Step 8: Check reciprocal liftover results and final cleanup----
species_consensus_categories = list()
for (sp in species){
  other_species = setdiff(species, sp)
  for (sp_other in other_species){
    filename = paste0(folder,sp,'_from_',sp_other, '_rd2_con_peaks')
    print(paste0('Processing ', sp, ' from ', sp_other,' peaks:'))
    result = concatenateAndSummarizeLiftovers(species_consensus_categories1[[sp]][[sp_other]],filename, peak_half_width, indel_high)
    reciprocal_peaks = result[[1]]; genome(reciprocal_peaks) = genome_names[[sp]]
    species_consensus_categories[[sp]][[sp_other]] <- checkReciprocalLiftovers(species_spec_peaks[[sp]],reciprocal_peaks,species_consensus_categories1[[sp]][[sp_other]], indel_values, summit_half_width_rd2, sp_other,genome_names)
  }
}

#Add peaks in consensus categories to consensus peaks
species_allcon_peaks1 = species_rd1_peaks
for (sp in species){
  species_allcon_peaks1[[sp]]$category = 'round1_indel_low' #existing peaks are round 1
}
for (sp in species){
  other_species = setdiff(species,sp)
  max_peak_name = species_allcon_peaks1[[sp]]$name[length(species_allcon_peaks1[[sp]])]; max_peak_num =  as.numeric(gsub("[^0-9]", "", max_peak_name)) + 1
  species_allcon_peaks1 = addMultiCorrespondingSpeciesPeaks(species_consensus_categories[[sp]],species_spec_peaks[[sp]], max_peak_num, species_allcon_peaks1, sp, other_species)
}
saveRDS(species_allcon_peaks1,paste0(folder,'species_allcon_peaks1.rds'))

#Exclude round 2 peaks that overlap round 1 peaks by more than overlap_tol or have lifted to blacklisted regions or chr_to_remove chromosomes
species_allcon_peaks2 <- excludeOverlappingAndBlacklistPeaks(species_allcon_peaks1, species_rd1_peaks,species_blacklists, species, remove_non_standard_chr, chr_to_remove, overlap_tol)
saveRDS(species_allcon_peaks2, paste0(folder,'species_allcon_peaks2.rds'))
#Exclude round 2 peaks that overlap other round 2 peaks more than tolerance from all species
species_allcon_peaks3 = species_allcon_peaks2
for (sp in species){#Iterate through species, keeping the first instance and removing all others from all species (first round should remove most from all species but there may be a few left due to imperfect liftovers so continue through all species)
  cleaned_peaks <- excludeSelfOverlaps(species_allcon_peaks3[[sp]],overlap_tol, ignore_strand = TRUE, keep = 'first')
  peaks_to_remove = setDiffForPeakNames(species_allcon_peaks3[[sp]],cleaned_peaks)
  peak_names_to_remove = peaks_to_remove$name
  for (sp in species){
    species_allcon_peaks3[[sp]] = species_allcon_peaks3[[sp]][species_allcon_peaks3[[sp]]$name %!in% peak_names_to_remove]
  }
}
saveRDS(species_allcon_peaks3, paste0(folder,'species_allcon_peaks3.rds'))
species_allcon_peaks = species_allcon_peaks3
for (sp in species){ #Remove metadata columns that were for internal use only
  species_allcon_peaks[[sp]]$category_orig <- NULL
  species_allcon_peaks[[sp]]$pass <- NULL
  species_allcon_peaks[[sp]]$summit_diff <- NULL
}
saveRDS(species_allcon_peaks, paste0(output_folder,'species_allcon_peaks.rds'))

#Add other categories to species_only_peaks, making sure peak names are still unique compared to those in species_allcon_peaks
species_only_peaks1 <- vector(mode = 'list', length = length(species)); names(species_only_peaks1) <- species
max_peak_name = species_allcon_peaks1[[sp]]$name[length(species_allcon_peaks1[[sp]])]; max_peak_num =  as.numeric(gsub("[^0-9]", "", max_peak_name)) + 1
for (i in 1:length(species)){
  sp = species[i]; other_species = setdiff(species,sp)
  species_only_peaks1[[sp]]$category = 'round1_indel_low' 
  species_only_peaks1[[sp]] = addMultiSpeciesOnlyPeaks(species_specific_categories[[sp]], species_spec_peaks[[sp]],max_peak_num, sp, other_species)
  #Make sure that peaks have unique names across species
  max_peak_name = species_only_peaks1[[sp]]$name[length(species_only_peaks1[[sp]])]; max_peak_num =  as.numeric(gsub("[^0-9]", "", max_peak_name)) + 1
}
saveRDS(species_only_peaks1,paste0(folder,'species_only_peaks1.rds'))

#Exclude peaks in species-only sets that overlap consensus sets by more than overlap tolerance (can happen when only part of a peak lifts over successfully)
species_only_peaks = species_only_peaks1
for (sp in species){
  overlapping_peak_names = findPeaksThatOverlapTooMuch(species_only_peaks[[sp]], species_allcon_peaks[[sp]], overlap_tol)
  species_only_peaks[[sp]] = species_only_peaks1[[sp]][species_only_peaks1[[sp]]$name %!in% overlapping_peak_names]
}
saveRDS(species_only_peaks,paste0(output_folder,'species_only_peaks.rds'))

#Export bed files
for (sp in species){
  export.bed(species_allcon_peaks[[sp]],paste0(output_folder,sp,'_allcon_peaks.bed'))
  export.bed(species_only_peaks[[sp]],paste0(output_folder,sp,'_only_peaks.bed'))
}

#Step 8: Final summaries----
#Summary of where peaks were lost in round 1
summary_steps1and2 <- summarizeConsensusPeaksPart1(species,peaks_orig,peaks,species_on_ancestral_summits1,species_on_ancestral_summits2,
                                         species_on_ancestral_peaks1,species_on_ancestral_peaks2)
print(format(summary_steps1and2,scientific=FALSE,digits = 2,nsmall = 2))
saveRDS(summary_steps1and2,paste0(output_folder,'summary_steps1and2.rds'))

summary_steps4and5 <- summarizeMultiConsensusPeaksPart2(species,consensus_peaks, consensus_df, species_con_summits1,species_con_summits2,
                                              species_rd1_peaks1,species_rd1_peaks2,species_rd1_peaks)
print(format(summary_steps4and5,scientific=FALSE,digits = 2,nsmall = 2))
saveRDS(summary_steps4and5,paste0(output_folder,'summary_steps4and5.rds'))

#Final summary
num_species <- length(species)
type_of_peak <- c('Consensus peaks rd 1', 'Consensus peaks rd 2','Total consensus')
num_con_peaks_rd1 <- length(species_rd1_peaks[[1]])
num_con_peaks_rd2 <- length(species_allcon_peaks[[1]]) - length(species_rd1_peaks[[1]])
num_con_peaks_total = length(species_allcon_peaks[[1]])
number_of_peaks <- c(num_con_peaks_rd1, num_con_peaks_rd2, num_con_peaks_total)
for (sp in species) {
  species_label <- paste(sp, "species-only peaks")
  type_of_peak <- c(type_of_peak, species_label)
  num_species_only_peaks <- length(species_only_peaks[[sp]])
  number_of_peaks <- c(number_of_peaks, num_species_only_peaks)
}
summary <- data.frame(type_of_peak, number_of_peaks)
saveRDS(summary,paste0(output_folder,'summary_final.rds'))
print(summary)

#Final checks (optional)----
#Check for no duplicate peaks
for (sp in species){
  all_peaks = c(species_allcon_peaks[[sp]],species_only_peaks[[sp]])
  mcols(all_peaks)$range_id <- paste0(seqnames(all_peaks), "-", start(all_peaks), "-", end(all_peaks))
  range_id_table <- table(mcols(all_peaks)$range_id)
  duplicate_range_ids <- names(range_id_table[range_id_table > 1])
  all_duplicate_indices <- which(mcols(all_peaks)$range_id %in% duplicate_range_ids)
  if (length(all_duplicate_indices)==0){print(paste0('Check successful - no duplicate ranges in final ',sp,' peak sets.'))
  } else {stop('Check failed - there are duplicate peaks in the final ',sp,' peak sets.')}
}

#Check no blacklisted or chr_to_remove peaks in final lists
no_blacklisted_peaks = TRUE
for (sp in species){
  result1 = findSingleOverlaps(species_allcon_peaks[[sp]],species_blacklists[[sp]])
  if (length(result1 > 0)){
    no_blacklisted_peaks = FALSE
  }
  result2 = findSingleOverlaps(species_only_peaks[[sp]],species_blacklists[[sp]])
  if (length(result2 > 0)){
    no_blacklisted_peaks = FALSE
  }
}
if (no_blacklisted_peaks) {
  print('Check successful - no blacklisted peaks are included in final lists')
} else {stop('Check failed - some blacklisted peaks are included in final lists')}

# Check if all species in species_allcon_peaks have the same number of peaks
all_equal_length <- all(sapply(2:length(species_allcon_peaks), function(i) length(species_allcon_peaks[[1]]) == length(species_allcon_peaks[[i]])))
if (all_equal_length) {
  print('Check successful - same number of peaks in species_allcon_peaks for all species')
} else {stop ('Check failed - species have different numbers of consensus peaks in species_allcon_peaks')}

# Check if all species in species_allcon_peaks have identical peak names
all_peaks_names <- lapply(species_allcon_peaks, function(x) x$name)
identical_names <- Reduce(function(x, y) x && identical(all_peaks_names[[1]], y), all_peaks_names, init = TRUE)
if (identical_names) {
  print('Check successful - all species have same names for consensus peaks')
} else {stop('Check failed - species have different names for consensus peaks')}

# Check if peak width is exactly set width for all peaks and only peaks
width_check_all <- all(sapply(species_allcon_peaks, function(peaks) max(width(peaks)) == 2*peak_half_width+1 && min(width(peaks)) == 2*peak_half_width+1))
width_check_only <- all(sapply(species_only_peaks, function(peaks) max(width(peaks)) == 2*peak_half_width+1 && min(width(peaks)) == 2*peak_half_width+1))
if (width_check_all & width_check_only) {
  print('Check successful - all peaks are the correct length')
} else {stop('Check failed - some peaks are not the same length')}

# Check if all peak names are unique for species only peaks
unique_names_all <- all(sapply(species_allcon_peaks, function(peaks) length(unique(peaks$name)) == length(peaks$name)))
unique_names_only <- all(sapply(species_only_peaks, function(peaks) length(unique(peaks$name)) == length(peaks$name)))
if (unique_names_all & unique_names_only) {
  print('Check successful - species_only_peaks have unique names')
} else {stop('Check failed - species_only peaks do not have unique names')}

# Check for no overlapping peaks between different species_only_peaks
overlap_free <- TRUE
for (i in 1:(length(species_only_peaks) - 1)) {
  for (j in (i + 1):length(species_only_peaks)) {
    if (!identical(intersect(species_only_peaks[[i]]$name, species_only_peaks[[j]]$name), character(0))) {
      overlap_free <- FALSE
      break
    }
  }
  if (!overlap_free) break
}
if (overlap_free) {
  print('Check successful - no redundant species-only peaks across species')
} else {stop('Check failed - some species_only_peaks do not have unique names')}

# Check for no overlapping peaks between species-only and consensus peaks
con_overlap_free <- TRUE
for (i in 1:length(species_only_peaks)) {
  if (!identical(intersect(species_only_peaks[[i]]$name, species_rd1_peaks[[i]]$name), character(0))) {
    con_overlap_free <- FALSE
    break
  }
}
if (con_overlap_free) {
  print('Check successful - no overlap between species-only and consensus peaks')
} else {stop('Check failed - some peaks are in consensus and species-only lists')}


#Percent of original peaks preserved
sp_percent_retained <- vector(mode = 'list', length = length(species)); names(sp_percent_retained) <- species
for (sp in species){
  original_peaks <- mcols(peaks[[sp]])$name
  species_peaks_orig_names <- mcols(species_allcon_peaks[[sp]])$orig_name
  split_orig_names <- strsplit(species_peaks_orig_names, ";")
  peaks_in_species_allcon <- unique(unlist(lapply(split_orig_names, function(x) {
    grep(paste0("^", prefixes[[sp]],"_peak_"), x, value = TRUE)
  })))
  peaks_in_species_only <- species_only_peaks[[sp]]$orig_name
  num_original_peaks <- length(original_peaks)
  num_retained_peaks_con <- sum(original_peaks %in% peaks_in_species_allcon)
  num_retained_peaks_only <- sum(original_peaks %in% peaks_in_species_only)
  sp_percent_retained[[sp]] <- ((num_retained_peaks_con + num_retained_peaks_only) / num_original_peaks) * 100
}
print('Percent of original peaks retained for each species: ') #note that the numbers are expected to be much lower for species beyond the first 2 since this doesn't account for consensus peaks that overlapped those
print(sp_percent_retained)