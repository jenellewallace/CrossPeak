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

#Setup blacklists----
human_blacklist_filename = 'human_blacklist.bed'
chimp_blacklist_filename = 'chimp_blacklist.bed'
rhesus_blacklist_filename = 'rhesus_blacklist.bed'
species_blacklists <- vector(mode = 'list',length = length(species)); names(species_blacklists) <- species
for (sp in species){
  if (sp == 'human'){
    species_blacklists[[sp]] = import.bed(paste0(folder,human_blacklist_filename))
  } else if (sp == 'chimp'){
    species_blacklists[[sp]] = import.bed(paste0(folder,chimp_blacklist_filename))
  } else if (sp == 'rhesus'){
    species_blacklists[[sp]] = import.bed(paste0(folder,rhesus_blacklist_filename))
  } else {stop('Unsupported species - no blacklist available. To use CrossPeak with this species, you will need to create your own blacklist and update the section Setup blacklists in CrossPeak_Step1.R')
  }
  genome(species_blacklists[[sp]]) = genome_names[[sp]]
} 
saveRDS(species_blacklists,paste0(folder,'species_blacklists.rds'))

#Step 1: Read peak file, make summits, and export summit bed files for halLiftover----
peaks_orig <- vector(mode = 'list', length = length(species)); names(peaks_orig) <- species
peaks <- vector(mode = 'list', length = length(species)); names(peaks) <- species
summits = vector(mode = 'list', length = length(species)); names(summits) <- species
for (sp in species){
  filename = paste0(sp,'_peaks')
  if (peaks_format == 'bed') {peaks_orig[[sp]] = import(paste0(folder,filename,'.bed'))
  } else if (peaks_format == 'rds'){peaks_orig[[sp]] = readRDS(paste0(folder,filename,'.rds'))
  } else {stop('Species peak sets must be in .rds or .bed format')}
  genome(peaks_orig[[sp]]) = genome_names[[sp]]
  #Make sure that peaks are all the same width (and width should be odd so that summit is perfect centered)
  if(max(width(peaks_orig[[sp]]))!=min(width(peaks_orig[[sp]]))){
    stop('This pipeline only works for fixed width peaks with summits centered and half width given by peak_half_width.')
  }
  if(!is_odd(max(width(peaks_orig[[sp]])))){
    warning('Peaks must have odd widths so summits can be perfectly centered. Subtracting one base pair from starts to ensure odd width.')
    current_starts = start(peaks_orig[[sp]])
    new_starts = current_starts - 1
    peaks_orig[[sp]]  <- GRanges(seqnames = seqnames(peaks_orig[[sp]]), 
                                 ranges = IRanges(start = new_starts, end = end(peaks_orig[[sp]])), 
                                 strand = strand(peaks_orig[[sp]]))
    genome(peaks_orig[[sp]]) = genome_names[[sp]]
  }
  if (keep_orig_peak_names & is.null(peaks_orig[[sp]]$name)){
    warning('The original peaks file must include a column called "name" if you want to retain the original peak names. No name column detected, defaulting to naming peaks with species prefixes.')
  }
  if (!keep_orig_peak_names | is.null(peaks_orig[[sp]]$name)){
    peaks[[sp]] <- labelPeaksWithPrefixes(peaks_orig[[sp]],sp, prefixes, keep.extra.columns = FALSE)
  } else {
    peaks[[sp]] <- peaks_orig[[sp]]
  }
  peaks[[sp]] <- suppressWarnings(excludeBlacklistPeaks(peaks[[sp]],species_blacklists[[sp]]))
  if (remove_non_standard_chr){
    non_std_chr = findNonStandardChromosomes(peaks[[sp]])
    peaks[[sp]] = dropSeqlevels(peaks[[sp]],non_std_chr, pruning.mode = 'coarse')
  }
  peaks[[sp]] = dropSeqlevels(peaks[[sp]],chr_to_remove, pruning.mode = 'coarse')
  strand(peaks[[sp]]) <- '+' #strand information not captured in peak calling but needed to keep track of strand flipping between original and ancestral genomes
  summits[[sp]] <- makeSummitsFromPeaks(peaks[[sp]], peak_half_width,summit_half_width)
  export.bed(summits[[sp]],paste0(folder,sp,'_summits.bed'))
}
saveRDS(peaks_orig, paste0(folder,'peaks_orig.rds'))
saveRDS(peaks, paste0(folder,'peaks.rds'))
saveRDS(summits, paste0(folder,'summits.rds'))