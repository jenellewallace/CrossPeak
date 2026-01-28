#Custom functions for working with GRanges objects
#should contain only the functions used in multi species pipeline version

#Adding back round 2 potentially species-specific peaks that succeeded in individual liftovers with their respective categories to the consensus set
addMultiCorrespondingSpeciesPeaks <- function(sp_con_categories,species_spec_peaks, max_peak_num, species_all_peaks, sp, other_species){
  #Figure out which peaks passed reciprocal liftover and are consensus for this species
  initial_species <- other_species[1]
  peak_names_to_add <- sp_con_categories[[initial_species]]$name[sp_con_categories[[initial_species]]$pass==TRUE]
  for (j in 1:length(other_species)) { #Find which peaks passed for all species
    sp_other = other_species[[j]]
    peak_names_to_add <- intersect(peak_names_to_add, sp_con_categories[[sp_other]]$name[sp_con_categories[[sp_other]]$pass==TRUE])
  }
  peakNumbers <- as.numeric(sub(".*_", "", peak_names_to_add)); orderIndex <- order(peakNumbers)
  orig_names <- peak_names_to_add[orderIndex]
  peaks_to_add <- vector(mode = "list", length = length(other_species)); names(peaks_to_add) = other_species
  for (j in 1:length(other_species)) { #Pull out coordinates for peaks that will be added to consensus list
    sp_other = other_species[[j]]
    peaks_to_add[[sp_other]] <- sp_con_categories[[sp_other]][sp_con_categories[[sp_other]]$name %in% orig_names]
  }
  #Check that they are properly sorted - this is crucial since they will be renamed and must correspond across species
  #all species should be in the same order because they were sorted by original peak name going in
  isSorted1 = is_sorted_by_peak_name(species_spec_peaks); 
  isSorted2 = TRUE
  for (sp_other in other_species){
    if (!is_sorted_by_peak_name(peaks_to_add[[sp_other]])){isSorted2 = FALSE}
  }
  if (!(isSorted1 & isSorted2)){stop('At least one of the input granges is not sorted.')}
  categories = '' #List of categories for sp lifted to each other species (in original order of al species)
  for (sp_other in other_species){
    categories = paste0(categories, ';', peaks_to_add[[sp_other]]$category)
  }
  new_names = paste0('Peak_',seq(from = max_peak_num, to = max_peak_num+length(orig_names)-1, by = 1))
  sp_coords = species_spec_peaks[species_spec_peaks$name %in% orig_names]
  sp_coords$orig_name = orig_names
  sp_coords$name = new_names
  names(sp_coords) = sp_coords$name
  sp_coords$category = categories
  sp_coords$category_orig = categories
  sp_coords$pass = TRUE #just for housekeeping, automatically passes since peak originated in this species
  species_all_peaks[[sp]] = combineGrs(species_all_peaks[[sp]], sp_coords)
  for (sp_other in other_species){
    sp_other_coords = peaks_to_add[[sp_other]]
    sp_other_coords$orig_name = orig_names
    sp_other_coords$name = new_names
    names(sp_other_coords) = sp_other_coords$name
    sp_other_coords$category = categories
    species_all_peaks[[sp_other]] = combineGrs(species_all_peaks[[sp_other]], sp_other_coords)
  }
  return(species_all_peaks)
}

#Making species specific sets
addMultiSpeciesOnlyPeaks <- function(species_specific_categories, species_spec_peaks,max_peak_num, sp, other_species){
  #Figure out which peaks are shared between categories in all other species
  initial_species <- other_species[1]
  peak_names_to_add <- species_specific_categories[[initial_species]]$name
  for (j in 1:length(other_species)) { #Find which peaks are in this category for all species
    sp_other = other_species[[j]]
    peak_names_to_add <- intersect(peak_names_to_add, species_specific_categories[[sp_other]]$name)
  }
  peakNumbers <- as.numeric(sub(".*_", "", peak_names_to_add)); orderIndex <- order(peakNumbers)
  orig_names <- peak_names_to_add[orderIndex]
  peaks_to_add <- vector(mode = "list", length = length(other_species)); names(peaks_to_add) = other_species
  for (j in 1:length(other_species)) { #Pull out coordinates for peaks that will be added to consensus list
    sp_other = other_species[[j]]
    peaks_to_add[[sp_other]] <- species_specific_categories[[sp_other]][species_specific_categories[[sp_other]]$name %in% orig_names]
  }
  categories = '' #List of categories for sp lifted to each other species (in original order)
  for (sp_other in other_species){
    categories = paste0(categories, ';', peaks_to_add[[sp_other]]$category)
  }
  new_names = paste0('Peak_',seq(from = max_peak_num, to = max_peak_num+length(orig_names)-1, by = 1))
  sp_coords = species_spec_peaks[species_spec_peaks$name %in% orig_names]
  sp_coords$orig_name = orig_names
  sp_coords$name = new_names
  names(sp_coords) = sp_coords$name
  sp_coords$category = categories
  species_only_peaks = sp_coords
  return(species_only_peaks)
}

#Assigns peaks to annotation categories based on which category the peak overlaps with the most (with intergenic being the default category for all
#ranges not specifically annotated). Ensure that duplicate ranges are removed from annotations object before using this function. Annotations object
#should have a column category with values belonging to category_types
assign_peaks_to_annotation_categories <- function(peaks, annotations, peak_width, category_types){
  # Find overlaps between peaks and annotations
  overlaps <- suppressWarnings(findOverlaps(peaks, annotations, ignore.strand = TRUE))
  # Extract indices of overlapping peaks and annotations
  peak_hits <- queryHits(overlaps)
  annotation_hits <- subjectHits(overlaps)
  # Get overlapping ranges
  overlapping_ranges <- pintersect(peaks[peak_hits], annotations[annotation_hits])
  # Compute overlap widths
  overlap_widths <- as.numeric(width(overlapping_ranges))
  # Extract categories as character vector
  categories <- as.character(mcols(annotations)$category)
  overlap_categories <- categories[annotation_hits]
  # Create data frame with standard column types
  df <- data.frame(
    peak_idx = as.integer(peak_hits),
    category = overlap_categories,
    overlap_width = overlap_widths,
    stringsAsFactors = FALSE
  )
  
  # Sum overlap widths within each category for each peak
  df_sum <- df %>%
    group_by(peak_idx, category) %>%
    summarise(total_overlap = sum(overlap_width), .groups = 'drop') %>%
    mutate(
      peak_idx = as.integer(peak_idx),
      category = as.character(category),
      total_overlap = as.numeric(total_overlap)
    )
  # Calculate total annotated overlap for each peak
  total_annotated_overlap <- df_sum %>%
    group_by(peak_idx) %>%
    summarise(total_annotated_overlap = sum(total_overlap), .groups = 'drop')
  # Compute intergenic overlap
  all_peak_idxs <- seq_along(peaks)
  df_total_overlap <- data.frame(peak_idx = as.integer(all_peak_idxs), stringsAsFactors = FALSE) %>%
    left_join(total_annotated_overlap, by = 'peak_idx') %>%
    mutate(
      total_annotated_overlap = ifelse(is.na(total_annotated_overlap), 0, total_annotated_overlap),
      intergenic_overlap = peak_width - total_annotated_overlap
    )
  # Include intergenic category
  df_intergenic <- df_total_overlap %>%
    mutate(
      category = 'intergenic',
      total_overlap = intergenic_overlap
    ) %>%
    dplyr::select(peak_idx, category, total_overlap)
  # Ensure consistent column types before combining
  df_sum <- df_sum %>%
    mutate(
      peak_idx = as.integer(peak_idx),
      category = as.character(category),
      total_overlap = as.numeric(total_overlap)
    )
  df_intergenic <- df_intergenic %>%
    mutate(
      peak_idx = as.integer(peak_idx),
      category = as.character(category),
      total_overlap = as.numeric(total_overlap)
    )
  # Combine all categories
  df_all_categories <- bind_rows(df_sum, df_intergenic)
  # For each peak, find the category with the maximum overlap
  df_max_overlap <- df_all_categories %>%
    group_by(peak_idx) %>%
    dplyr::filter(total_overlap == max(total_overlap)) %>%
    summarise(
      category = first(category),
      total_overlap = first(total_overlap),
      .groups = 'drop'
    )
  # Apply overlap threshold
  df_max_overlap <- df_max_overlap %>%
    mutate(
      category = ifelse(total_overlap > peak_width/2, category, 'intergenic')
    )
  # Summarize counts
  category_counts <- df_max_overlap %>%
    group_by(category) %>%
    summarise(count = n(), .groups = 'drop') %>%
    right_join(
      data.frame(category = category_types, stringsAsFactors = FALSE),
      by = 'category'
    ) %>%
    mutate(count = ifelse(is.na(count), 0, count))
  if (sum(category_counts$count) != length(peaks)){stop('Not all peaks were annotated. Check that annotations object includes category column.')}
  return(category_counts)
}

#Adds unique peaks for each species after the first two to theconsensus list and excludes peaks that overlap more than tolerance
addUniquePeaks <- function(existingPeaks, newPeaks, overlap_tol, next_species) {
  # Find overlaps with existing peaks
  overlaps = findOverlaps(newPeaks, existingPeaks, minoverlap = overlap_tol, ignore.strand = TRUE)
  
  # Identify non-overlapping new peaks
  nonOverlappingIndices = setdiff(seq_along(newPeaks), queryHits(overlaps))
  nonOverlappingNewPeaks = newPeaks[nonOverlappingIndices]
  nonOverlappingNewPeaks$orig_name = nonOverlappingNewPeaks$name
  nonOverlappingNewPeaks$orig_species = next_species
  
  # Combine existing peaks with non-overlapping new peaks
  updatedPeaks = c(existingPeaks, nonOverlappingNewPeaks)
  return(updatedPeaks)
}

#Adds metadata about liftover failures
categorizeConsensusPeaksFailures <- function(consensus_peaks,peaks_list, sp, failure_code){
  if (length(peaks_list) > 0) {
    indices <- match(peaks_list, names(consensus_peaks))
    indices <- na.omit(indices) 
    # Update the liftback_result with vectorized operations
    consensus_peaks$liftback_result[indices] <- ifelse(
      nzchar(consensus_peaks$liftback_result[indices]),
      paste0(consensus_peaks$liftback_result[indices], ';', failure_code),
      failure_code
    )
    # Update the liftback_species similarly
    consensus_peaks$liftback_species[indices] <- ifelse(
      nzchar(consensus_peaks$liftback_species[indices]),
      paste0(consensus_peaks$liftback_species[indices], ';', sp),
      sp
    )
  } else {
    warning(paste0('No ', failure_code, ' peaks for species ', sp))
  }
  return(consensus_peaks)
}

#Checks whether round 2 peaks lift back reciprocally to peak of origin with different tolerances depending on category
#Also changes category if liftover succeeds but indel tolerance is higher than it should be for that category
checkReciprocalLiftovers <- function(spec_peaks,reciprocal_peaks,con_categories, indel_values, summit_half_width, sp_other, genome_names){
  #Setup
  tolerance_values = (indel_values + 2 * summit_half_width) / 2 - summit_half_width #How much a summit would move under each indel condition
  get_new_category <- function(pass_current, pass_any, category, summit_diff) {
    if (pass_current) {
      return(category)
    } else if (!pass_current && pass_any) {
      # Explicitly check the current category and determine the next category based on summit_diff
      if (category == "round2_indel_low" && summit_diff < tolerance_values["round2_indel_med"]) {
        return("round2_indel_med")
      } else if (category == "round2_indel_low" && summit_diff < tolerance_values["round2_indel_high"]) {
        return("round2_indel_high")
      } else if (category == "round2_indel_med" && summit_diff < tolerance_values["round2_indel_high"]) {
        return("round2_indel_high")
      }
    }
    return(category)  # Default to the original category if no transitions apply
  }
  spec_peaks_df <- as.data.frame(spec_peaks)
  reciprocal_peaks_df <- as.data.frame(reciprocal_peaks)
  con_categories_df <- as.data.frame(con_categories)
  #Find summits
  spec_peaks_df <- spec_peaks_df %>%
    dplyr::mutate(summit = (start + end) / 2)
  reciprocal_peaks_df <- reciprocal_peaks_df %>%
    dplyr::mutate(summit = (start + end) / 2)
  #Join orig and reciprocal data frames and add info
  peaks_check <- inner_join(spec_peaks_df,reciprocal_peaks_df, by = 'name')
  cat_indices <- match(peaks_check$name,con_categories_df$name)
  peaks_check$category = con_categories_df$category[cat_indices]
  peaks_check <- peaks_check %>%
    dplyr::mutate(seqnames.x = as.character(seqnames.x),
                  seqnames.y = as.character(seqnames.y)) %>%
    dplyr::mutate(tolerance = (indel_values[category] + 2 * summit_half_width) / 2 - summit_half_width) %>%
    dplyr::mutate(summit_diff = abs(summit.x - summit.y)) %>%
    dplyr::mutate(pass_current = ifelse(seqnames.x == seqnames.y & summit_diff < tolerance, TRUE, FALSE)) %>%
    dplyr::mutate(pass_any = ifelse(seqnames.x == seqnames.y & summit_diff < max(tolerance), TRUE, FALSE)) %>%
    dplyr::mutate(category_new = mapply(get_new_category, pass_current, pass_any, category, summit_diff))
  #Add info back to consensus_categories 
  match_indices = match(con_categories_df$name,peaks_check$name)
  con_categories_df$category_new = peaks_check$category_new[match_indices]
  con_categories_df$pass = peaks_check$pass_any[match_indices]
  con_categories_df$pass[is.na(con_categories_df$pass)] = FALSE  #set false for peaks that didn't match (ie failed or multiple reciprocal lift)
  con_categories_df$summit_diff = peaks_check$summit_diff[match_indices]
  con_categories_new <- GRanges(
    seqnames = con_categories_df$seqnames,
    ranges = IRanges(start = con_categories_df$start, end = con_categories_df$end),
    strand = con_categories_df$strand,
    name = con_categories_df$name,
    category_orig = con_categories_df$category,
    category = con_categories_df$category_new,
    pass = con_categories_df$pass,
    summit_diff = con_categories_df$summit_diff
  )
  genome(con_categories_new) = genome_names[[sp_other]]
  return(con_categories_new)
}


#Combine granges objects, checking whether genomes are the same
#(default behavior only throws error if they are different but not if one is 
#unspecified whereas here this function throws an error in both of these cases)
combineGrs <- function(...){
  grs_list <- list(...)
  
  # Extract genomes from each GenomicRanges object and flatten them into a single vector
  genomes_vector <- unlist(lapply(grs_list, function(gr) unique(genome(gr))))
  
  # Check if all genomic sequences refer to the same genome
  if (length(unique(genomes_vector)) != 1) {
    stop('Genomic ranges objects cannot be combined because not all refer to the same genome.')
  } else {
    # Use do.call with c to combine all GenomicRanges objects in the list
    result <- do.call(c, grs_list)
    return(result)
  }
}

#Combines all regions with the same name on the same chromosome into a single region and returns granges with
#combined regions and data frame with the names of the regions that failed due to mapping to regions on different chromosomes, 
#Input must be a granges object with a name column
combineLiftoverFragments <- function(gr) {
  # Convert GenomicRanges to a data frame for easier manipulation
  gr_df <- as.data.frame(gr)
  # Group by 'name' and summarize each group
  combined <- gr_df %>%
    group_by(.data$name) %>%
    summarize(
      seqnames = if(n_distinct(.data$seqnames) > 1) NA_character_ else .data$seqnames[[1]],
      start = min(.data$start),
      end = max(.data$end),
      strand = if(n_distinct(.data$strand) > 1) NA_character_ else .data$strand[[1]],
      .groups = 'drop'
    )
  # Split the combined data frame into successes and failures - failure is when peaks of the same name mapped to different chromosomes and/or strands
  successes <- combined %>% dplyr::filter(!is.na(.data$seqnames) & !is.na(.data$strand))
  failures <- combined %>% dplyr::filter(is.na(.data$seqnames) | is.na(.data$strand))
  # Convert successes back to GRanges
  peaks_primate_lo_combined <- makeGRangesFromDataFrame(successes, keep.extra.columns = TRUE)
  peaks_primate_lo_combined_sorted <- sortGrByName(peaks_primate_lo_combined)
  # Return combined GRanges and failures
  return(list(peaks_primate_lo_combined_sorted, failures))
}

#Check if final peaks are in the same place as original peaks (after lifting to ancestral and back), considering movement of the summit for merged peaks
compareMultiConsensusAndOriginalSummits <- function(orig_summits, merge_metadata, species_con_summits, indel_tol, summit_half_width, prefix) {
  # Ensure the data is in the correct format for dplyr operations
  merge_metadata <- as.data.frame(merge_metadata)
  orig_summits <- as.data.frame(orig_summits)
  species_con_summits <- as.data.frame(species_con_summits)
  #Extract only peaks for species of interest
  species_con_summits <- species_con_summits %>%
    mutate(peak_of_interest = str_extract(.data$orig_name, paste0(prefix, "_peak_\\d+")))
  filtered_species_con_summits <- species_con_summits %>%
    filter(!is.na(.data$peak_of_interest))
  names(filtered_species_con_summits)[names(filtered_species_con_summits) == "peak_of_interest"] <- paste0(prefix, "_peak_of_origin")
  
  # Create a dataframe from merge_metadata
  df <- merge_metadata
  orig_summits_adjusted <- orig_summits %>% 
    dplyr::rename(original_name = .data$name)
  # Perform the join
  joined_summits <- df %>%
    left_join(orig_summits_adjusted, by = c("original_peak_name" = "original_name")) %>%
    mutate(
      # Calculate the original summit adjustment
      original_summit_adj = ifelse(!is.na(.data$start), as.numeric(.data$start) + summit_half_width + as.numeric(.data$summit_adj), NA),
      # Calculate the consensus summit position
      con_summit_summit = species_con_summits$start[match(.data$consensus_name, species_con_summits$name)] + summit_half_width,
      # Calculate the difference between the consensus summit and the original summit adjustment
      summit_diffs = .data$con_summit_summit - .data$original_summit_adj,
      # Calculate a dynamic tolerance based on the adjustment
      current_tols = ifelse(abs(as.numeric(.data$summit_adj))/10 < 1, indel_tol, round(abs(as.numeric(.data$summit_adj)/10) * indel_tol))
    ) 
  
  problem_summits <- joined_summits %>%
    select(.data$consensus_name, .data$original_peak_name,.data$summit_adj, .data$summit_diffs, .data$current_tols)
  
  # Filtering problem summits
  problem_summits_final <- problem_summits %>%
    filter(!is.na(.data$summit_diffs), !is.na(.data$current_tols), abs(.data$summit_diffs) > .data$current_tols)
  
  print('Problem summits identified.')
  return(problem_summits_final)
}

#Concatenates liftover results of the same names and prints a summary
concatenateAndSummarizeLiftovers <- function(human_summits, filename, summit_half_width, indel_tol){
  lo <- import(paste0(filename,'_lo.bed'))
  result <- combineLiftoverFragments(lo)
  human_peaks_lo_all = result[[1]]; human_peaks_lo_multiple = result[[2]]
  names(human_peaks_lo_all) = human_peaks_lo_all$name
  summary_result = summarizeLiftoverResults(human_summits,human_peaks_lo_all,human_peaks_lo_multiple, summit_half_width)
  summary = summary_result$summary
  saveRDS(summary, paste0(filename,'_summary.rds'))
  human_peaks_lo <- excludeDifferentWidthsSummits(human_peaks_lo_all, indel_tol, summit_half_width)
  names(human_peaks_lo) = human_peaks_lo$name  
  return(list(all_lifted_peaks = human_peaks_lo_all,lifted_peaks_exclude_with_indel_tol = human_peaks_lo, 
              failed_lifted_peaks = human_peaks_lo_multiple, summary_result = summary_result))
}

#Excludes peaks overlapping blacklist region from a single species peak set
excludeBlacklistPeaks <- function(peaks,current_blacklist){
  blacklist_peaks = findSingleOverlaps(peaks,current_blacklist)
  peaks = peaks[!peaks$name %in% blacklist_peaks$name,]
  return(peaks)
}

#Separates out summits that had a different width than the original after liftover (they can be included later depending on indel_tol)
excludeDifferentWidthsSummits <- function(gr, indel_tol, half_width){
  widths = width(gr)
  valid_indices <- which(widths >= 2*half_width+1 - indel_tol & widths <= 2*half_width+1 +indel_tol)
  new_gr = gr[valid_indices,]
  seqlevels(new_gr) <- sort(seqlevels(new_gr))
  new_gr <- sort(new_gr)
}

#Add columns to consensus_peaks describing peaks that failed to lift back to individual species genomes and exclude these peaks from all species consensus peaks
excludeFailedProblemBlacklistPeaks<- function(consensus_peaks,con_peaks,problem_peaks, blacklist, species, remove_non_standard_chr, chr_to_remove){
  #List the names of peaks that did not lift back successfully for each category and note in consensus_peaks
  failed_peaks_list <- vector(mode = 'list', length = length(species)); names(failed_peaks_list) = species
  problem_peaks_list <- vector(mode = 'list', length = length(species)); names(problem_peaks_list) = species
  blacklist_peaks_list <- vector(mode = 'list', length = length(species)); names(blacklist_peaks_list) = species
  chr_remove_peaks_list <- vector(mode = 'list', length = length(species)); names(chr_remove_peaks_list) = species
  consensus_peaks$liftback_result = character(length(consensus_peaks))
  consensus_peaks$liftback_species = character(length(consensus_peaks))
  con_peaks2 = con_peaks
  for (sp in species){
    matched_indices <- match(mcols(consensus_peaks)$name, mcols(con_peaks[[sp]])$name)
    #Failed peaks
    has_no_match <- is.na(matched_indices)
    failed_peaks_names = names(consensus_peaks)[has_no_match]
    failed_peaks_list[[sp]] = failed_peaks_names
    consensus_peaks <- categorizeConsensusPeaksFailures(consensus_peaks,failed_peaks_list[[sp]], sp, 'Failed')
    #Problem peaks
    problem_peaks_list[[sp]] = problem_peaks[[sp]]$consensus_name
    consensus_peaks <- categorizeConsensusPeaksFailures(consensus_peaks,problem_peaks_list[[sp]], sp, 'Problem')
    #Blacklist peaks
    result = findSingleOverlaps(con_peaks[[sp]],blacklist[[sp]])
    blacklist_peaks_list[[sp]] = result$name
    consensus_peaks <- categorizeConsensusPeaksFailures(consensus_peaks,blacklist_peaks_list[[sp]], sp, 'Blacklist')
    #Peaks on non-standard and excluded chromosomes
    if (remove_non_standard_chr){
      non_std_chr = findNonStandardChromosomes(con_peaks2[[sp]])
      con_peaks2[[sp]] = dropSeqlevels(con_peaks2[[sp]],non_std_chr, pruning.mode = 'coarse')
    }
    con_peaks2[[sp]] = dropSeqlevels(con_peaks2[[sp]],chr_to_remove, pruning.mode = 'coarse')
    peaks_removed = setdiff(con_peaks[[sp]]$name,con_peaks2[[sp]]$name)
    chr_remove_peaks_list[[sp]] = peaks_removed
    consensus_peaks <- categorizeConsensusPeaksFailures(consensus_peaks,chr_remove_peaks_list[[sp]], sp, 'Chr_removed')
  }
  failed_peaks = unique(unlist(failed_peaks_list))
  problem_peaks = unique(unlist(problem_peaks_list))
  blacklist_peaks = unique(unlist(blacklist_peaks_list))
  chr_remove_peaks = unique(unlist(chr_remove_peaks_list))
  
  #Exclude peaks from all species
  exclude_peaks = c(failed_peaks,problem_peaks,blacklist_peaks, chr_remove_peaks)
  con_peaks_final <- con_peaks
  for (sp in species){
    con_peaks_final[[sp]] = con_peaks[[sp]][!con_peaks[[sp]]$name %in% exclude_peaks,]
  }
  print('Failed, problem, and blacklisted peaks removed from all species consensus peaks.')
  if (remove_non_standard_chr){'Peaks on non-standard chromosomes removed from all species consensus peaks.'}
  if (is.null(chr_to_remove)){'Peaks on chr_to_remove removed from all species consensus peaks.'}
  result = list(consensus_peaks,con_peaks_final)
}

#Exclude round 2 peaks that overlap more than tolerance with round 1 peaks as well as blacklisted or chr_to_remove peaks
excludeOverlappingAndBlacklistPeaks <- function(species_all_peaks1, species_con_peaks,
                                                blacklist, species, remove_non_standard_chr, chr_to_remove, overlap_tol){
  overlapping_peak_names <- vector(mode = 'list', length = length(species)); names(overlapping_peak_names) <- species
  blacklist_peak_names  <- vector(mode = 'list', length = length(species)); names(blacklist_peak_names) <- species
  chr_to_remove_peak_names  <- vector(mode = 'list', length = length(species)); names(chr_to_remove_peak_names) <- species
  species_all_peaks_test = species_all_peaks1
  for (sp in species){
    #Overlapping peaks - exclude round2 overlapping round1
    new_peaks = setDiffForPeakNames(species_all_peaks1[[sp]],species_con_peaks[[sp]])
    overlapping_peak_names[[sp]] = findPeaksThatOverlapTooMuch(new_peaks, species_con_peaks[[sp]], overlap_tol)
    #Blacklist peaks
    result = findSingleOverlaps(species_all_peaks1[[sp]],blacklist[[sp]])
    blacklist_peak_names[[sp]] = result$name
    #Peaks on non-standard and excluded chromosomes
    if (remove_non_standard_chr){
      non_std_chr = findNonStandardChromosomes(species_all_peaks_test[[sp]])
      species_all_peaks_test[[sp]] = dropSeqlevels(species_all_peaks_test[[sp]],non_std_chr, pruning.mode = 'coarse')
    }
    species_all_peaks_test[[sp]] = dropSeqlevels(species_all_peaks_test[[sp]],chr_to_remove, pruning.mode = 'coarse')
    peaks_removed = setdiff(species_all_peaks1[[sp]]$name,species_all_peaks_test[[sp]]$name)
    chr_to_remove_peak_names[[sp]] = peaks_removed  
  }
  #Exclude peaks
  all_peaks_to_remove = c(unique(unlist(overlapping_peak_names)),unique(unlist(blacklist_peak_names)), unique(unlist(chr_to_remove_peak_names)))
  species_all_peaks2 = species_all_peaks1
  for (sp in species){
    species_all_peaks2[[sp]] = species_all_peaks2[[sp]][!(species_all_peaks2[[sp]]$name %in% all_peaks_to_remove)] 
  }
  return(species_all_peaks2)
}

#Excludes peaks that overlap within the same granges object (e.g. after liftover) by more than tol with option to exclude all or keep the first instance
excludeSelfOverlaps <- function(gr, tol, ignore_strand, keep = 'none') {
  # Find overlaps within the same GRanges object
  overlaps <- findOverlaps(gr, gr, minoverlap = tol, ignore.strand = ignore_strand)
  valid_overlaps <- overlaps[queryHits(overlaps) != subjectHits(overlaps)]
  
  # Get indices
  overlapping_indices <- unique(c(queryHits(valid_overlaps), subjectHits(valid_overlaps)))
  non_overlapping_indices <- setdiff(seq_along(gr), subjectHits(valid_overlaps))
  
  # Subset the GRanges object to keep only non-overlapping ranges
  if (keep == 'none') {
    gr_filtered <- gr[non_overlapping_indices]
  }
  
  if (keep == 'first') {
    # Create a data frame of overlapping pairs
    overlap_pairs <- data.frame(
      query = queryHits(valid_overlaps),
      subject = subjectHits(valid_overlaps)
    )
    
    # Ensure the first instance is the one with the smaller start coordinate
    first_instance <- overlap_pairs %>%
      mutate(min_start = pmin(start(gr[query]), start(gr[subject])),
             pair_id = pmin(query, subject) * 1e10 + pmax(query, subject)) %>%
      group_by(pair_id) %>%
      slice_min(min_start, with_ties = FALSE) %>%
      pull(query)
    
    # Combine first instance indices with non-overlapping indices
    gr_filtered <- gr[sort(unique(c(first_instance, non_overlapping_indices)))]
  }
  num_removed = length(gr) - length(gr_filtered)
  print(paste0(num_removed, ' peaks were excluded because they overlapped by more than ', tol, ' bp.'))
  return(gr_filtered)
}

#Makes a new granges object that expands a granges object by the defined amount
expandRange = function(x, upstream, downstream){
  strand_is_minus = strand(x) == "-"
  on_plus = which(!strand_is_minus)
  on_minus = which(strand_is_minus)
  start(x)[on_plus] = start(x)[on_plus] - upstream
  start(x)[on_minus] = start(x)[on_minus] - downstream
  end(x)[on_plus] = end(x)[on_plus] + downstream
  end(x)[on_minus] = end(x)[on_minus] + upstream
  x }

#Filters a granges object to retain only peaks with only user-selected categories. User must select which species categories should be used for filtering
gr_filter <- filter_for_high_confidence <- function(gr, categories, species_index){
  #Split category column by semiocolon
  get_entry_by_index <- function(category, col_index) {
    entries <- unlist(strsplit(category, ";"))
    # Remove empty entries and select the one at the desired index
    entries <- entries[entries != ""]
    if (length(entries) >= col_index) {
      return(entries[col_index])
    } else {
      return(NA)  # Return NA if the index is out of bounds
    }
  }
  col_index = species_index - 1
  d1 = as.data.frame(gr)
  d2 = d1 %>%   # Apply the get_entry_by_index function to extract the desired category for each row
    mutate(sp_category = sapply(category, get_entry_by_index, col_index = col_index))
  # Filter based on user-selected categories
  d3 = d2[d2$sp_category %in% categories,]
  # Convert back to GRanges object
  gr_filter = makeGRangesFromDataFrame(d3, keep.extra.columns = TRUE)
  return(gr_filter)
}

#Finds non-standard chromosomes in primate genomes - may need to modify for other species
findNonStandardChromosomes <- function(gr) {
  # Extract the sequence names from the GRanges object
  seq_names <- seqnames(gr)
  
  # Define a pattern that matches non-standard chromosomes
  # This pattern matches strings that do not start with 'chr'
  # or have an underscore followed by additional characters after 'chrN'
  pattern <- "^((?!chr).+|chr[0-9XY]+_.+|chrUn.+)$"
  
  # Find and return sequence names that match the pattern
  non_standard_chrs <- seq_names[grepl(pattern, seq_names, perl = TRUE)]
  non_standard_chrs = as.character(unique(non_standard_chrs))
  return(non_standard_chrs)
}


#Finds all overlaps between gr1 and gr2 and returns lengths
findOverlapLengths <- function(gr1,gr2, ignore_strand, min_overlap = 0){
  commonSeqlevels <- intersect(seqlevels(gr1), seqlevels(gr2))
  gr1 <- gr1[seqnames(gr1) %in% commonSeqlevels]
  gr2 <- gr2[seqnames(gr2) %in% commonSeqlevels]
  overlaps <- findOverlaps(gr1, gr2, ignore.strand = ignore_strand, minoverlap = min_overlap)
  overlapLengths <- pmin(end(gr1[queryHits(overlaps)]), end(gr2[subjectHits(overlaps)])) -
    pmax(start(gr1[queryHits(overlaps)]), start(gr2[subjectHits(overlaps)])) + 1
  return(overlapLengths)
}

#Prints a summary of overlaps between two granges objects and returns the overlaps
findOverlapDetails <- function(gr1, gr2, ignore_strand,species_names, print_result = TRUE) {
  commonSeqlevels <- intersect(seqlevels(gr1), seqlevels(gr2))
  #gr1 <- gr1[seqnames(gr1) %in% commonSeqlevels]
  #gr2 <- gr2[seqnames(gr2) %in% commonSeqlevels]
  
  overlaps <- findOverlaps(gr1, gr2, ignore.strand = ignore_strand)
  overlapLengths <- pmin(end(gr1[queryHits(overlaps)]), end(gr2[subjectHits(overlaps)])) -
    pmax(start(gr1[queryHits(overlaps)]), start(gr2[subjectHits(overlaps)])) + 1
  
  # Count of unique overlapping ranges in gr1 and gr2
  numOverlaps <- list(length(unique(queryHits(overlaps))),length(unique(subjectHits(overlaps))))
  
  # Count of non-overlapping ranges in gr1 and gr2
  numNonOverlapsGr1 <- length(gr1) - numOverlaps[[1]]
  numNonOverlapsGr2 <- length(gr2) - numOverlaps[[2]]
  numNonOverlaps <- list(numNonOverlapsGr1,numNonOverlapsGr2)
  
  labels = c('Non-overlapping peaks','Peaks with overlaps','Total')
  df = data.frame(labels, species1 = c(c(numNonOverlaps[[1]],numOverlaps[[1]],length(gr1))),
                  species2 = c(numNonOverlaps[[2]],numOverlaps[[2]],length(gr2)))
  names(df)[2:3] <- species_names
  if(print_result){print(df)}
  
  return(list(overlapLengths = overlapLengths, 
              numOverlaps = numOverlaps, 
              numNonOverlaps = numNonOverlaps))
}

#Finds peaks in gr1 that overlap more than tolerance with those in gr2 and return list of their names
findPeaksThatOverlapTooMuch <- function(gr1,gr2, overlap_tol){
  overlaps <- findOverlaps(gr1, gr2, ignore.strand = TRUE) 
  overlapLengths <- findOverlapLengths(gr1, gr2, TRUE)
  # Create a data frame to map overlaps
  overlapMap <- data.frame(gr1Index = queryHits(overlaps),
                           gr2Index = subjectHits(overlaps),
                           overlapLen = overlapLengths)
  # Filter overlaps based on the overlap tolerance
  validOverlaps <- overlapMap[overlapMap$overlapLen >= overlap_tol, ]
  overlap_names = gr1$name[validOverlaps$gr1Index]
  return(overlap_names)
}

#Find overlaps within a single granges object
findSelfOverlaps <- function(gr){
      # Find overlaps within the same GRanges object
      overlaps <- findOverlaps(gr, gr, ignore.strand = TRUE)
      valid_overlaps <- overlaps[queryHits(overlaps) != subjectHits(overlaps)]
      # Get indices of ranges that are not involved in valid overlaps
      overlapLengths <- pmin(end(gr[queryHits(valid_overlaps)]), end(gr[subjectHits(valid_overlaps)])) -
        pmax(start(gr[queryHits(valid_overlaps)]), start(gr[subjectHits(valid_overlaps)])) + 1
      return(overlapLengths)
}

#Find overlaps from one gr (ex. peaks) to another (ex. a blacklist) and return a gr containing the rows of the overlapping peaks from the first list
findSingleOverlaps <- function(gr1, gr2) {
  overlaps <- findOverlaps(gr1, gr2)
  overlap_regions = gr1[queryHits(overlaps)]
  return(overlap_regions)
}

#Takes two granges objects with a names column, finds which regions in each overlap 
#2 or more regions in the #other and then outputs two data frames showing the 
#correspondences
formatOverlaps <- function(query, subject, overlaps) {
  queryNames <- mcols(query)$name[queryHits(overlaps)]
  subjectNames <- mcols(subject)$name[subjectHits(overlaps)]
  overlapsList <- split(subjectNames, queryNames)
  # Filter out entries with less than two overlaps
  overlapsList <- overlapsList[sapply(overlapsList, length) >= 2]
  overlapsList <- lapply(overlapsList, function(x) paste(x, collapse=", "))
  overlapsList
}

#For a list of peak names in one species, gets the cooresponding coordinates from species_all_peaks and species_only_peaks
getCoordsForPeakNames <- function(peak_names,sp,species_all_peaks, species_only_peaks){
  all_peaks = combineGrs(species_all_peaks[[sp]], species_only_peaks[[sp]])
  available_peaks = names(all_peaks)
  peak_names = gsub("-", "_", peak_names);
  valid_peak_names = peak_names[peak_names %in% available_peaks]
  coords = all_peaks[valid_peak_names]
  return(coords)
}

#Given two granges objects with names columns describing matching regions, return a data frame with the widths of the 
#matching regions (missing regions will be discarded)
getMatchingRangeWidths <- function(gr1,gr2){
  widths_human <- width(gr1)
  if (is.null(gr1$name)){stop('Missing name column')}
  names(widths_human) <- gr1$name
  widths_primate_lo <- width(gr2)
  if (is.null(gr2$name)){stop('Missing name column')}
  names(widths_primate_lo) <- gr2$name
  matched_ids <- intersect(gr1$name, gr2$name)
  matched_widths_human <- widths_human[matched_ids]
  matched_widths_primate_lo <- widths_primate_lo[matched_ids]
  df <- data.frame(
    WidthsGr1 = matched_widths_human,
    WidthsGr2 = matched_widths_primate_lo
  )
  return(df)
}

#Gets putative species-specific peaks (peaks that didn't lift to ancestral or that failed to lift back to other species), will be tested further in round 2
getMultiSpeciesSpecificPeaks <- function(species, species_on_primate_peaks2, peaks, consensus_peaks, consensus_df){
  species_orig <- vector(mode = 'list', length = length(species)); names(species_orig) <- species
  merged <- vector(mode = 'list', length = length(species)); names(merged) <- species
  species_spec_peaks1 <- vector(mode = 'list', length = length(species)); names(species_spec_peaks1) <- species
  fail_codes = c('Failed','Problem','Blacklist','Chr_removed')
  for (sp in species){
    #Part1
    if (!is.null(species_on_primate_peaks2[[sp]]$name)){
      part1_names = setdiff(peaks[[sp]]$name,species_on_primate_peaks2[[sp]]$name) 
      part1 = peaks[[sp]][peaks[[sp]]$name %in% part1_names] #peaks originating from species i that failed to lift to ancestral
    } else {stop('Names are null.')}
    
    #Part 2
    fail_indices = which(consensus_df$liftback_result %in% fail_codes & consensus_df$orig_species == sp)
    # fail_names = unique(consensus_df$name[fail_indices])
    # part2[[sp]] = consensus_peaks[fail_names] #peaks originating from species i that failed to lift back from ancestral to any species
    
    #Get original coords
    part1_orig_names = part1$name
    part2_orig_names = consensus_df$orig_name[fail_indices]
    orig_names = c(part1_orig_names,part2_orig_names)
    matching_indices <- match(peaks[[sp]]$name, orig_names)
    species_spec_peaks1[[sp]] =peaks[[sp]][!is.na(matching_indices), ]
  }
  #Sort by name
  species_spec_peaks <- vector(mode = 'list', length = length(species)); names(species_spec_peaks) <- species
  for (sp in species){
    species_spec_peaks[[sp]] = sortGrByName(species_spec_peaks1[[sp]])
    strand(species_spec_peaks[[sp]]) <-  '*' #not considering strand info for quantification
  }
  print('Potential species-specific peaks identified for follow-up.')
  return(species_spec_peaks)
}

#Finds peaks overlapping promoters
get_promoter_peaks <- function(peak_granges, annotation, upstream = 1000, downstream = 100){
  tss <- GetTSSPositions(annotation)
  promoters = expandRange(tss, upstream = upstream, downstream = downstream)
  # Find overlaps between peaks and promoters
  overlaps <- findOverlaps(peak_granges, promoters)
  overlapping_peaks <- queryHits(overlaps)
  promoter_peaks_ranges <- peak_granges[overlapping_peaks]
  return(promoter_peaks_ranges)
}

#Makes a granges object from astring with a dash separating parts
makeGrangesFromString <- function(ranges_string) {
  parts <- strsplit(ranges_string, "-")[[1]]
  chromosome <- parts[1]
  start <- as.numeric(parts[2])
  end <- as.numeric(parts[3])
  GRanges(seqnames = chromosome, ranges = IRanges(start = start, end = end))
}

#Gets the names of peaks from a Signac object within a specified range
getPeaksFromPlot <- function(obj,assay,extend_up,extend_down, gene = NULL, range_string = NULL){
  annotation = Annotation(obj)
  DefaultAssay(obj) = assay
  peaks_celltypes = granges(obj)
  if (!is.null(gene)){
    gene_range <- subset(annotation, gene_name == gene)
    extended_range <- GRanges(
      seqnames = seqnames(gene_range),
      ranges = IRanges(
        start = start(gene_range) - extend_up,
        end = end(gene_range) + extend_down
      ),
      strand = strand(gene_range))
  } else if (!is.null(range)){
    range = makeGrangesFromString(range_string)
    extended_range <- GRanges(
      seqnames = seqnames(range),
      ranges = IRanges(
        start = start(range) - extend_up,
        end = end(range) + extend_down
      ),
      strand = strand(range))
  }
  #seqlevels(extended_range) <- paste0("chr", seqlevels(extended_range)) # Rename the sequence levels to have the "chr" prefix
  overlapping_peaks <- subsetByOverlaps(peaks_celltypes, extended_range)
  return(overlapping_peaks)
}

#Gets the names of peaks from a Granges object within a specified range
getPeaksinRange <- function(peaks_celltypes,annotation,extend_up,extend_down, gene = NULL, range_string = NULL){
  if (!is.null(gene)){
    gene_range <- subset(annotation, gene_name == gene)
    extended_range <- GRanges(
      seqnames = seqnames(gene_range),
      ranges = IRanges(
        start = start(gene_range) - extend_up,
        end = end(gene_range) + extend_down
      ),
      strand = strand(gene_range))
  } else if (!is.null(range)){
    range = makeGrangesFromString(range_string)
    extended_range <- GRanges(
      seqnames = seqnames(range),
      ranges = IRanges(
        start = start(range) - extend_up,
        end = end(range) + extend_down
      ),
      strand = strand(range))
  }
  #seqlevels(extended_range) <- paste0("chr", seqlevels(extended_range)) # Rename the sequence levels to have the "chr" prefix
  overlapping_peaks <- subsetByOverlaps(peaks_celltypes, extended_range)
  return(overlapping_peaks)
}

#Categorizes round2 peaks based on separate summit and peak liftovers and indel tolerances
getRound2PeakCategories <- function(species,species_spec_peaks, species_spec_summits, indel_low, indel_med, indel_high, 
                                    genome_names,lift_folder, peak_half_width, summit_half_width){
  num_species <- length(species); num_other_species <- num_species-1
  empty_list <- vector("list", num_species * num_other_species)
  #Peaks where the summit lifted with varying indel tolerances - will be added to consensus peaks
  species_perfect_peaks <- vector("list", num_species); names(species_perfect_peaks) = species
  species_indel_low_peaks<- vector("list", num_species); names(species_indel_low_peaks) = species
  species_indel_med_peaks<- vector("list", num_species); names(species_indel_med_peaks) = species
  species_indel_high_peaks<- vector("list", num_species); names(species_indel_high_peaks) = species
  #Peaks that failed to lift entirely
  species_failed_peaks <- vector("list", num_species); names(species_failed_peaks) = species
  #Peaks where summit failed but whole peak succeeded with width less than peak width
  species_summit_deletion_peaks <- vector("list", num_species); names(species_summit_deletion_peaks) = species
  #Peaks where summit lifted but width was so large that summit would move outside of existing peak - set to 0 for other species in count matrix
  species_summit_with_large_indel_peaks <- vector("list", num_species); names(species_summit_with_large_indel_peaks) = species
  #Peaks where summit failed and width was so large that summit would move outside of existing peak - set to 0 for other species in count matrix
  species_summit_failed_large_indel_peaks <- vector("list", num_species); names(species_summit_failed_large_indel_peaks) = species
  #Consolidated categories
  species_consensus_categories <- vector("list", num_species); names(species_consensus_categories) = species
  species_specific_categories <- vector("list", num_species); names(species_specific_categories) = species
  
  for (sp in species){
    other_species = setdiff(species,sp)
    species_perfect_peaks[[sp]] <- vector("list", num_other_species); names(species_perfect_peaks[[sp]]) = other_species
    species_indel_low_peaks[[sp]] <- vector("list", num_other_species); names(species_indel_low_peaks[[sp]]) = other_species
    species_indel_med_peaks[[sp]] <- vector("list", num_other_species); names(species_indel_med_peaks[[sp]]) = other_species
    species_indel_high_peaks[[sp]] <- vector("list", num_other_species); names(species_indel_high_peaks[[sp]]) = other_species
    species_failed_peaks[[sp]] <- vector("list", num_other_species); names(species_failed_peaks[[sp]]) = other_species
    species_summit_deletion_peaks[[sp]] <- vector("list", num_other_species); names(species_summit_deletion_peaks[[sp]]) = other_species
    species_summit_with_large_indel_peaks[[sp]] <- vector("list", num_other_species); names(species_summit_with_large_indel_peaks[[sp]]) = other_species
    species_summit_failed_large_indel_peaks[[sp]] <- vector("list", num_other_species); names(species_summit_failed_large_indel_peaks[[sp]]) = other_species
    species_consensus_categories[[sp]] <- vector("list", num_other_species); names(species_consensus_categories[[sp]]) = other_species
    species_specific_categories[[sp]] <- vector("list", num_other_species); names(species_specific_categories[[sp]]) = other_species
    for (sp_other in other_species){
      num_invalid = 0 #counter for peaks that were excluded due to ranges exceeding limits
      filename = paste0(lift_folder,sp,'_on_',sp_other,'_spec_summits')
      print(paste0('Processing ', sp, ' on ', sp_other,' peaks:'))
      lo <- import(paste0(filename,'_lo.bed'))
      result <- combineLiftoverFragments(lo)
      summits_lo_all = result[[1]]; summits_lo_multiple = result[[2]]
      names(summits_lo_all) = summits_lo_all$name
      genome(summits_lo_all) = genome_names[[sp_other]]
      result = summarizeLiftoverResults(species_spec_summits[[sp]],summits_lo_all,summits_lo_multiple, summit_half_width)
      saveRDS(result$summary, paste0(filename,'_summary.rds'))
      summit_failed_names <- result$names_failed
      #Perfect summits - add to consensus peaks as is (low indel_tol)
      species_perfect_summits = summits_lo_all[summits_lo_all$name %in% result$names_same]
      summits_lo1 = setDiffForPeakNames(summits_lo_all,species_perfect_summits)
      #Summits within indel_tol 5 - add to consensus peaks as is (low_indel_tol)
      species_indel_low_summits <- excludeDifferentWidthsSummits(summits_lo1, indel_low, summit_half_width)
      summits_lo2 = setDiffForPeakNames(summits_lo1,species_indel_low_summits)
      #Summits within indel_tol 50 - add to consensus peaks with med_indel_tol
      species_indel_med_summits <- excludeDifferentWidthsSummits(summits_lo2, indel_med, summit_half_width)
      summits_lo3 = setDiffForPeakNames(summits_lo2,species_indel_med_summits)
      #Summits within indel_tol 250 - add to consensus peaks with high_indel_tol
      species_indel_high_summits <- excludeDifferentWidthsSummits(summits_lo3, indel_high, summit_half_width)
      summits_lo4 = setDiffForPeakNames(summits_lo3,species_indel_high_summits) 
      species_summit_with_large_indel_peaks[[sp]][[sp_other]] = species_spec_peaks[[sp]][species_spec_peaks[[sp]]$name %in% summits_lo4$name] #get original peaks for that species
      species_summit_with_large_indel_peaks[[sp]][[sp_other]] = sortGrByName(species_summit_with_large_indel_peaks[[sp]][[sp_other]]); species_summit_with_large_indel_peaks[[sp]][[sp_other]]$category = 'round2_summit_with_large_indel'
      #Check that all peaks are accounted for
      num_peaks = length(summits_lo4) + length(species_indel_high_summits) +length(species_indel_med_summits) + length(species_indel_low_summits) +length(species_perfect_summits)
      if (num_peaks != length(summits_lo_all)){stop('Number of successfully lifted potential species-specific summits does not add up.')}
      #Rebuild peaks around summits
      species_perfect_peaks[[sp]][[sp_other]] <- makePeaksFromSummitRanges(species_perfect_summits,peak_half_width)
      species_perfect_peaks[[sp]][[sp_other]] = sortGrByName(species_perfect_peaks[[sp]][[sp_other]]); species_perfect_peaks[[sp]][[sp_other]]$category = 'round2_indel_low'
      num_invalid = num_invalid + length(species_perfect_summits) - length(species_perfect_peaks[[sp]][[sp_other]])
      species_indel_low_peaks[[sp]][[sp_other]] <- makePeaksFromSummitRanges(species_indel_low_summits,peak_half_width)
      species_indel_low_peaks[[sp]][[sp_other]] = sortGrByName(species_indel_low_peaks[[sp]][[sp_other]]); species_indel_low_peaks[[sp]][[sp_other]]$category = 'round2_indel_low'
      num_invalid = num_invalid + length(species_indel_low_summits) - length(species_indel_low_peaks[[sp]][[sp_other]])
      species_indel_med_peaks[[sp]][[sp_other]] <- makePeaksFromSummitRanges(species_indel_med_summits,peak_half_width)
      species_indel_med_peaks[[sp]][[sp_other]] = sortGrByName(species_indel_med_peaks[[sp]][[sp_other]]); species_indel_med_peaks[[sp]][[sp_other]]$category = 'round2_indel_med'
      num_invalid = num_invalid + length(species_indel_med_summits) - length(species_indel_med_peaks[[sp]][[sp_other]])
      species_indel_high_peaks[[sp]][[sp_other]] <- makePeaksFromSummitRanges(species_indel_high_summits,peak_half_width)
      species_indel_high_peaks[[sp]][[sp_other]] <- sortGrByName(species_indel_high_peaks[[sp]][[sp_other]]); species_indel_high_peaks[[sp]][[sp_other]]$category = 'round2_indel_high'
      num_invalid = num_invalid + length(species_indel_high_summits) - length(species_indel_high_peaks[[sp]][[sp_other]])
      #Get peaks that failed to lift over entirely
      filename = paste0(lift_folder,sp,'_on_',sp_other,'_spec_peaks')
      lo <- import(paste0(filename,'_lo.bed'))
      result <- combineLiftoverFragments(lo)
      peaks_lo_all = result[[1]]; peaks_lo_multiple = result[[2]]
      names(peaks_lo_all) = peaks_lo_all$name
      genome(peaks_lo_all) = genome_names[[sp_other]]
      result = summarizeLiftoverResults(species_spec_peaks[[sp]],peaks_lo_all,peaks_lo_multiple, peak_half_width)
      species_peaks_failed_summits = species_spec_peaks[[sp]][species_spec_peaks[[sp]]$name %in% result$names_failed]
      species_failed_peaks[[sp]][[sp_other]] <- makePeaksFromSummitRanges(species_peaks_failed_summits,peak_half_width)
      species_failed_peaks[[sp]][[sp_other]] <- sortGrByName(species_failed_peaks[[sp]][[sp_other]]); species_failed_peaks[[sp]][[sp_other]]$category = 'round2_peak_failed'
      num_invalid = num_invalid + length(species_peaks_failed_summits) - length(species_failed_peaks[[sp]][[sp_other]])
      saveRDS(result$summary, paste0(filename,'_summary.rds'))
      
      #Species summit deletion peaks
      #summit failed but whole peak succeeded
      summit_failed_peak_success = intersect(summit_failed_names,peaks_lo_all$name)
      species_summit_deletion_all = peaks_lo_all[peaks_lo_all$name %in% summit_failed_peak_success]
      species_summit_deletion_summits = species_summit_deletion_all[width(species_summit_deletion_all)<2*peak_half_width+1]
      species_summit_deletion_peaks[[sp]][[sp_other]] = makePeaksFromSummitRanges(species_summit_deletion_summits,peak_half_width)
      species_summit_deletion_peaks[[sp]][[sp_other]] <- sortGrByName(species_summit_deletion_peaks[[sp]][[sp_other]]); species_summit_deletion_peaks[[sp]][[sp_other]]$category = 'round2_summit_deletion'
      num_invalid = num_invalid + length(species_summit_deletion_summits) - length(species_summit_deletion_peaks[[sp]][[sp_other]])
      species_summit_failed_large_indel_summits_othersp = setDiffForPeakNames(species_summit_deletion_all,species_summit_deletion_summits) 
      species_summit_failed_large_indel_peaks[[sp]][[sp_other]] = species_spec_peaks[[sp]][species_spec_peaks[[sp]]$name %in% species_summit_failed_large_indel_summits_othersp$name] #get original peaks for that species
      species_summit_failed_large_indel_peaks[[sp]][[sp_other]]<- sortGrByName(species_summit_failed_large_indel_peaks[[sp]][[sp_other]]); species_summit_failed_large_indel_peaks[[sp]][[sp_other]]$category = 'round2_summit_failed_large_indel'
      #summit failed and whole peak was multi-mapping or summit and peak were both multimapping
      summit_failed_peak_multi = intersect(summit_failed_names, peaks_lo_multiple$name)
      summit_multi_peak_multi = intersect(summits_lo_multiple$name, peaks_lo_multiple$name)
      all_multi = c(summit_failed_peak_multi,summit_multi_peak_multi)
      total_multi = length(unique(all_multi))
      
      #Consolidate into consensus and species-specific categories
      species_consensus_categories[[sp]][[sp_other]] = combineGrs(species_perfect_peaks[[sp]][[sp_other]],species_indel_low_peaks[[sp]][[sp_other]],species_indel_med_peaks[[sp]][[sp_other]],
                                                                  species_indel_high_peaks[[sp]][[sp_other]],species_summit_deletion_peaks[[sp]][[sp_other]])
      species_consensus_categories[[sp]][[sp_other]] = sortGrByName(species_consensus_categories[[sp]][[sp_other]])
      species_specific_categories[[sp]][[sp_other]] = combineGrs(species_summit_with_large_indel_peaks[[sp]][[sp_other]], species_failed_peaks[[sp]][[sp_other]],
                                                                 species_summit_failed_large_indel_peaks[[sp]][[sp_other]])
      species_specific_categories[[sp]][[sp_other]] = sortGrByName(species_specific_categories[[sp]][[sp_other]])
      
      #Summarize
      category = c('Peaks with summits that lifted perfectly', 'Peaks with summits that lifted with low indel tolerance', 
                   'Peaks with summits that lifted with med indel tolerance','Peaks with summits that lifted with high indel tolerance',
                   'Peaks with summits that lifted but with indels larger than peak width', 'Peaks that failed to lift in entirety', 
                   'Peaks with summit deletions', 'Peaks with summit deletions and large insertions elsewhere', 
                   'Peaks with multi-mapping summits or failed summits and multi-mapping peaks', 'Total')
      values = c(length(species_perfect_peaks[[sp]][[sp_other]]) , length(species_indel_low_peaks[[sp]][[sp_other]]) , length(species_indel_med_peaks[[sp]][[sp_other]]) ,
                 length(species_indel_high_peaks[[sp]][[sp_other]]) , length(species_summit_with_large_indel_peaks[[sp]][[sp_other]]) , length(species_failed_peaks[[sp]][[sp_other]]) , 
                 length(species_summit_deletion_peaks[[sp]][[sp_other]]) , length(species_summit_failed_large_indel_peaks[[sp]][[sp_other]]), total_multi, length(species_spec_summits[[sp]]))
      summary = data.frame(category,values); print(summary)
      
      #Check that all peaks are unique and no overlap between categories
      num_total = sum(values[1:(length(values)-1)]) + num_invalid
      if (length(species_spec_summits[[sp]]) != num_total){stop('Number of species-specific summits across all categories (',num_total, ') does not add up to the total.')}
      overlap = intersect(species_consensus_categories[[sp]][[sp_other]]$name, species_specific_categories[[sp]][[sp_other]]$name)
      if (length(overlap)!=0){stop('Some peaks are in multiple categories.')}
      saveRDS(summary,paste0(lift_folder, sp, '_species_specific_peaks_summary.rds'))
    }
  }
  return(list(species_consensus_categories = species_consensus_categories, species_specific_categories = species_specific_categories))
}

'%!in%' <- function(x,y)!('%in%'(x,y))

#Initializes metadata during creation of round 1 consensus peak set
initializeMetadataDataFrames <- function(speciesList, speciesNames) {
  metadataList <- list()
  for (i in seq_along(speciesList)) {
    gr <- speciesList[[i]]
    df <- data.frame(
      consensus_name = character(length(gr)),
      original_peak_name = mcols(gr)$name,
      original_chr = seqnames(gr),
      original_start = start(gr),
      original_end = end(gr),
      strand_flip = strand(gr) == "-",
      summit_adj = numeric(length(gr)),  # Initially zero
      stringsAsFactors = FALSE
    )
    metadataList[[speciesNames[i]]] <- df
  }
  return(metadataList)
}

is_odd <- function(x) {
  if (x %% 2 != 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#Checks whether a granges object is sorted by the peak name
is_sorted_by_peak_name <- function(gr) {
  # Extract numeric part of the peak names
  peak_numbers <- as.numeric(gsub(".*_([0-9]+)$", "\\1", mcols(gr)$name))
  # Check if the numeric parts are in ascending order
  all(diff(peak_numbers) > 0)
}

#Creates a names column with named peaks for each species according to the defined prefixes
labelPeaksWithPrefixes <- function(summits, current_species,prefixes, keep.extra.columns = FALSE){
  current_prefix <- tryCatch({
    prefixes[[current_species]]
  }, error = function(e) {
    "Un"
  })
  if (!keep.extra.columns){
    mcols(summits) <- NULL
  }
  peak_names <- paste0(current_prefix, '_peak_', seq_along(summits))
  mcols(summits)$name <- peak_names
  return(summits)
}

#Makes fixed width peaks from granges objects If summit_col == NULL then it assumes summit is centered within summit range and extends 
#both sides to create a peak with given half width. Otherwise it takes the summit as the start position plus the value in summit_col
makePeaksFromSummitRanges <- function(gr, peak_half_width, summit_col = NULL) {
  if (is.null(summit_col)){
    summit_positions = mid(gr)
    starts <- summit_positions - peak_half_width
    ends <- summit_positions + peak_half_width
  } else {
    summit_positions = mcols(gr)[[summit_col]]
    starts <- start(gr) + summit_positions - peak_half_width
    ends <- start(gr) + summit_positions + peak_half_width
  }
  
  # Filter out peaks where the start would be negative
  valid_indices <- which(starts >= 1)
  
  if (length(valid_indices) > 0) {
    new_gr <- GRanges(seqnames = seqnames(gr)[valid_indices], 
                      ranges = IRanges(start = starts[valid_indices], end = ends[valid_indices]), 
                      strand = strand(gr)[valid_indices])
    # Assigning metadata columns
    mcols(new_gr) <- mcols(gr)[valid_indices,]
    if("X" %in% names(mcols(new_gr))) {
      names(mcols(new_gr))[names(mcols(new_gr)) == "X"] <- "name"
    }
    new_gr <- sort(new_gr)
    names(new_gr) <- new_gr$name
    genome(new_gr) = genome(gr)
  }
  num_invalid = length(which(starts < 1))
  print(paste0('Peaks made from summit ranges (', num_invalid,' summits had invalid starts.)'))
  return(new_gr)
}

#Makes summits with specified width based on peaks, assues summit is centered within peak
makeSummitsFromPeaks <- function(gr, peak_half_width, summit_half_width) {
  # Calculate the start and end positions of the summits
  summit_starts = start(gr) + peak_half_width - summit_half_width
  summit_ends = start(gr) + peak_half_width + summit_half_width
  
  # Create a new GRanges object for the summits
  summit_gr = GRanges(
    seqnames = seqnames(gr),
    ranges = IRanges(start = summit_starts, end = summit_ends),
    strand = strand(gr)
  )
  
  # Copy all metadata from the original gr to the new summit_gr
  mcols(summit_gr) = mcols(gr)
  genome(summit_gr) = genome(gr)
  names(summit_gr) = names(gr)
  return(summit_gr)
}

#Merges peaks across an arbitrary number of species. Merges the first two species based on two species merging rules
#and then for the subsequent species, discards any peaks that overlap the existing set by more than overlap_tol while
#adding in any unique peaks or peaks that overlap less than tolerance.
mergePeaksAcrossMultipleSpecies <- function(speciesList, speciesNames, prefixes, overlap_tol, half_width) {
  if (length(speciesList) < 2) stop("At least two species are required for merging.")
  
  # Initialize metadata for each species
  metadataList <- initializeMetadataDataFrames(speciesList, speciesNames)
  
  # Start with the peaks from the first species
  mergedPeaks <- speciesList[[1]]
  mergedPeaks$orig_name = speciesList[[1]]$name
  mergedPeaks$orig_species = speciesNames[[1]]
  
  # Merge the first two species
  first_two_species = speciesNames[1:2]
  mergedPeaks <- mergeTwoSpeciesPeaks(mergedPeaks, speciesList[[2]], first_two_species,overlap_tol, half_width)
  print(paste0('Successfully merged peaks from ',speciesNames[[1]], ' and ', speciesNames[[2]]))
  
  # Process each subsequent species
  if (length(speciesList)>2){
    for (i in 3:length(speciesList)) {
      current_species_gr <- speciesList[[i]]
      current_species_name <- speciesNames[[i]]
      mergedPeaks <- addUniquePeaks(mergedPeaks, current_species_gr, overlap_tol,current_species_name)
      print(paste0('Successfully merged peaks from ', current_species_name))
    }
  }
  
  # Sort and rename merged peaks
  mergedPeaks = sort(mergedPeaks)
  names(mergedPeaks) <- paste0("Peak_", seq_along(mergedPeaks))
  mcols(mergedPeaks)$name <- names(mergedPeaks)
  
  # Update metadata based on the final merged peaks
  metadataList <- updateMetadataPostMerge(mergedPeaks, metadataList, speciesList, speciesNames)
  
  return(list(mergedPeaks = mergedPeaks, metadata = metadataList))
}

#Merges a list of granges objects (e.g., from different cell types) into a single set of peaks (using macs2 p-values)
mergePeaksWithIterativeOverlap <- function(granges_list, num_cores, sig_col = 'neg_log10qvalue_summit') {
  # Apply removeOverlappingPeaks to each GRanges object
  list1 <- bplapply(granges_list, removeOverlappingPeaks, sig_col = sig_col, BPPARAM = MulticoreParam(workers = num_cores))
  print('Peaks merged within cell types.')
  
  # Merge all processed GRanges objects
  list2 = GRangesList(list1)
  merged_granges <- unlist(list2)
  rownames(merged_granges) <- NULL
  
  # Apply removeOverlappingPeaks to the merged GRanges object
  final_granges <- removeOverlappingPeaks(merged_granges, sig_col)
  
  # Sort the final GRanges object
  sorted_final_granges <- sort(final_granges)
  print('Peaks merged across cell types - iterative overlap merging complete.')
  
  return(sorted_final_granges)
}

#Merges peaks from the first two species
mergeTwoSpeciesPeaks <- function(gr1, gr2, two_species, overlap_tol, half_width) {
  overlaps = findOverlaps(gr1, gr2, minoverlap = overlap_tol, ignore.strand = TRUE)
  mergedRanges = gr1[setdiff(seq_along(gr1), queryHits(overlaps))]  # Start with non-overlapping from gr1
  
  if (length(overlaps) > 0) {
    idx1 = queryHits(overlaps)
    idx2 = subjectHits(overlaps)
    
    averageSummits = round(rowMeans(cbind(start(gr1[idx1]), end(gr2[idx2]))))
    newStarts = averageSummits - half_width
    newEnds = newStarts + 2 * half_width
    
    # Filtering valid new starts
    validIdx = newStarts >= 1
    newRanges = GRanges(
      seqnames = seqnames(gr1)[idx1][validIdx],
      ranges = IRanges(start = newStarts[validIdx], end = newEnds[validIdx]),
      strand = rep('*', sum(validIdx))
    )
    
    mcols(newRanges)$orig_name <- paste(gr1[idx1][validIdx]$name, gr2[idx2][validIdx]$name, sep=";")
    mcols(newRanges)$orig_species <- paste(two_species[[1]], two_species[[2]], sep=";")
    mergedRanges = c(mergedRanges, newRanges)
  }
  
  # Include non-overlapping ranges from gr2
  gr2_non_overlap_peaks = gr2[setdiff(seq_along(gr2), subjectHits(overlaps))]
  mcols(gr2_non_overlap_peaks)$orig_name <- gr2_non_overlap_peaks$name
  mcols(gr2_non_overlap_peaks)$orig_species <- two_species[[2]]
  mergedRanges = c(mergedRanges, gr2_non_overlap_peaks)
  
  return(mergedRanges)
}

#Generates and saves histograms of the overlap lengths between all pairs of species
plot_overlap_histograms <- function(species_on_ancestral_peaks, species, save_plots = TRUE, figure_folder){
  # Find overlap details for each combination of species
  results <- combn(seq_along(species_on_ancestral_peaks), 2, function(indices) {
    # Accessing the corresponding species names based on the current indices
    species1 <- species[indices[1]]
    species2 <- species[indices[2]]
    
    # Accessing the peak data from 'species_on_ancestral_peaks' by indices
    peaks1 <- species_on_ancestral_peaks[[indices[1]]]
    peaks2 <- species_on_ancestral_peaks[[indices[2]]]
    
    result <- findOverlapDetails(peaks1, peaks2, ignore_strand = TRUE, species = c(species1,species2), print_result = TRUE)
    list(
      overlap_lengths = result[[1]],
      numOverlaps = result[[2]],
      numNonOverlaps = result[[3]]
    )
  }, simplify = FALSE)
  
  # Naming the results
  names(results) <- combn(species, 2, function(sp) {
    paste("Comparison", sp[1], "vs", sp[2])
  }, simplify = TRUE)
  
  # Define histogram breaks
  breaks_values <- seq(from = 0, to = 520, by = 20)  # Adjust as necessary
  
  # Generate and save plots for each pair
  for (name in names(results)) {
    data <- data.frame(overlap_lengths = results[[name]]$overlap_lengths)
    p <- ggplot(data, aes(x = overlap_lengths)) +
      geom_histogram(breaks = breaks_values, color = "black", fill = "white") +
      theme_minimal() +
      labs(title = name,  # Use the name directly as the title
           x = "Length of overlap", y = "Number of cross-species peak pairs",
           caption = paste0('Total pairs of overlapping peaks: ', length(data$overlap_lengths),
                            ', Median: ', median(data$overlap_lengths)))
    
    if (save_plots){
      ggsave(filename = paste0(figure_folder, gsub(" ", "_", name), '_histogram_of_overlap_lengths.png'),
             plot = p, height = 4, width = 6, dpi = 300)
    }
  }
}

#Generates and saves histograms of the overlap lengths between all pairs of species
plot_overlap_self_histograms <- function(peaks, species, save_plots = TRUE, figure_folder){
  for (sp in species){
    gr = peaks[[sp]]
    # Find overlaps within the same GRanges object
    overlaps <- findOverlaps(gr, gr, ignore.strand = TRUE)
    valid_overlaps <- overlaps[queryHits(overlaps) != subjectHits(overlaps)]
    # Get indices of ranges that are not involved in valid overlaps
    overlapLengths <- pmin(end(gr[queryHits(valid_overlaps)]), end(gr[subjectHits(valid_overlaps)])) -
      pmax(start(gr[queryHits(valid_overlaps)]), start(gr[subjectHits(valid_overlaps)])) + 1
    
    # Generate and save plots for each species
    breaks_values <- seq(from = 0, to = 520, by = 20)  # Adjust as necessary
    data <- data.frame(overlap_lengths = overlapLengths)
    p <- ggplot(data, aes(x = overlap_lengths)) +
      geom_histogram(breaks = breaks_values, color = "black", fill = "white") +
      theme_minimal() +
      labs(title = str_to_title(sp),  # Use the name directly as the title
           x = "Length of overlap", y = "Number of overlapping peak pairs",
           caption = paste0('Total pairs of overlapping peaks: ', length(data$overlap_lengths),
                            ', Median: ', median(data$overlap_lengths)))
    
    if (save_plots){
      ggsave(filename = paste0(figure_folder, sp, '_histogram_of_overlap_lengths.png'),
             plot = p, height = 4, width = 6, dpi = 300)
    }
  }
}

#Removes peaks that overlap as part of mergePeaksWithIterativeOverlap
removeOverlappingPeaks <- function(granges_obj, sig_col) {
  # Ensure there are peaks to process
  if (length(granges_obj) == 0) {
    return(GRanges())
  }
  
  # Check if sig_col column exists and is numeric
  if (!sig_col %in% names(mcols(granges_obj))) {
    stop(paste0(sig_col," column is missing"))
  }
  if (!is.numeric(mcols(granges_obj)[[sig_col]])) {
    stop(paste0(sig_col, " column is not numeric"))
  }
  
  # Remove NA values if present
  valid_indices <- which(!is.na(mcols(granges_obj)[[sig_col]]))
  granges_obj <- granges_obj[valid_indices]
  
  # Initialize an empty GRanges object for retained peaks
  retained_peaks <- GRanges()
  
  while(length(granges_obj) > 0) {
    # Sort by significance, preserving original order in case of ties
    sorted_indices <- order(-mcols(granges_obj)[[sig_col]], seq_along(granges_obj))
    granges_obj <- granges_obj[sorted_indices]
    
    # Retain the most significant peak (the first one in case of ties)
    most_significant_peak <- granges_obj[1]
    retained_peaks <- c(retained_peaks, most_significant_peak)
    
    # Identify and remove overlapping peaks
    overlapping_indices <- subjectHits(findOverlaps(most_significant_peak, granges_obj))
    granges_obj <- granges_obj[-overlapping_indices]
  }
  retained_peaks
}

#Performs the setdiff function but using peak names
#Returns gr3 which includes only peaks that are in gr1 that are not in gr2 based on name column
setDiffForPeakNames <- function(gr1, gr2){
  names_to_exclude <- setdiff(gr1$name, gr2$name) 
  gr3 = gr1[gr1$name %in% names_to_exclude] #get object with perfect summits excluded for next step
  return(gr3)
}

#Sort granges object by name column
sortGrByName <- function(gr){
  peakNumbers <- as.numeric(sub(".*_", "", mcols(gr)$name))
  orderIndex <- order(peakNumbers)
  sortedGr <- gr[orderIndex]
  return(sortedGr)
}

# Within one granges object, finds overlapping peaks and adds half of range to each one
# Returns a granges object with same length as original and spanning same genomic regions but with nonoverlapping peaks
split_overlapping_peaks <- function(gr) {
  # Convert GRanges to data frame
  gr_df <- as.data.frame(gr)
  # Find overlaps
  overlaps <- findOverlaps(gr, ignore.strand = TRUE)
  overlap_pairs <- as.data.frame(overlaps)
  # Get unique overlapping pairs
  overlap_pairs <- overlap_pairs %>% 
    dplyr::filter(queryHits < subjectHits)
  if (nrow(overlap_pairs) == 0) return(gr)  # If no overlaps, return original GRanges
  # Extract overlapping start and end coordinates
  overlap_pairs <- overlap_pairs %>%
    dplyr::mutate(overlap_start = pmax(gr_df$start[queryHits], gr_df$start[subjectHits]),
           overlap_end = pmin(gr_df$end[queryHits], gr_df$end[subjectHits]),
           overlap_width = overlap_end - overlap_start + 1) %>%
    dplyr::filter(overlap_width > 1)  # Only keep meaningful overlaps
  # Calculate half overlap width
  overlap_pairs <- overlap_pairs %>%
    mutate(half_overlap_width = floor(overlap_width / 2))
  # Adjust the start and end positions in a vectorized way
  gr_df <- gr_df %>%
    mutate(end = ifelse(row_number() %in% overlap_pairs$queryHits,
                        end - overlap_pairs$half_overlap_width[match(row_number(), overlap_pairs$queryHits)], end),
           start = ifelse(row_number() %in% overlap_pairs$subjectHits,
                          start + overlap_pairs$half_overlap_width[match(row_number(), overlap_pairs$subjectHits)], start))
  # Convert the adjusted data frame back to GRanges
  gr_adjusted <- makeGRangesFromDataFrame(gr_df, keep.extra.columns = TRUE)
  return(gr_adjusted)
}

#Summarize where peaks were lost in round 1, part 1 - lifting to the ancestral genome
summarizeConsensusPeaksPart1 <- function(species,peaks_orig,peaks,species_on_primate_summits1,species_on_primate_summits2,
                                         species_on_primate_peaks1,species_on_primate_peaks2){
  #Part 1 - lifting to ancestral
  df_part1_sp = list()
  for (sp in species){
    df_step1 = length(peaks_orig[[sp]]) - length(peaks[[sp]]) #peaks overlapping blacklist or alternate chr
    df_stepA = length(peaks[[sp]]) - length(species_on_primate_summits1[[sp]]) #peaks that failed liftover (failed and multiple)
    df_step2a = length(species_on_primate_summits1[[sp]]) - length(species_on_primate_summits2[[sp]]) # peaks that failed liftover (different widths)
    df_step2b = length(species_on_primate_summits2[[sp]]) - length(species_on_primate_peaks1[[sp]]) # peaks that failed liftover (overlapping ends of chromosomes)
    df_step2c = length(species_on_primate_peaks1[[sp]]) - length(species_on_primate_peaks2[[sp]]) # peaks that failed liftover (overlapping)
    total_removed = length(peaks_orig[[sp]])- length(species_on_primate_peaks2[[sp]])
    #Check that total removed is equal to sum of removed at each step
    total_all_steps = df_step1+df_stepA+df_step2a+df_step2b+df_step2c
    if (total_all_steps == total_removed){print(paste0('Part 1 check successful - total peaks removed is equal to the sum of peaks removed from each step for species ',sp))
    } else {warning('Part 1 check failed - total peaks removed is NOT equal to the sum of peaks removed from each step for ', sp)}
    df_part1_sp[[sp]] = c(length(peaks_orig[[sp]]),df_step1,df_stepA,df_step2a, df_step2b,df_step2c,
                          total_removed,length(species_on_primate_peaks2[[sp]]),
                          length(species_on_primate_peaks2[[sp]])/length(peaks_orig[[sp]]))
  }
  df_part1 <- do.call(cbind, df_part1_sp)
  rownames(df_part1) = c('Total to start','Removed, Step 1','Removed, Step A', 'Removed, Step 2a',
                         'Removed, Step 2b','Removed, Step 2c','Total removed','Total remaining', 'Proportion remaining')
  return(df_part1)
}

#Summarizes the results of lifting over peaks from one species to another
#Returns a list containing the summary and named categories of peaks
summarizeLiftoverResults <- function(peaks, peaks_lo, peaks_lo_multiple, half_width, show = TRUE){
  total = length(peaks)
  lo = length(peaks_lo)
  df <- getMatchingRangeWidths(peaks,peaks_lo)
  df_same <- df[df$WidthsGr1 == df$WidthsGr2, ]
  df_dif <- df[df$WidthsGr1 != df$WidthsGr2, ]
  same = dim(df_same)[1]
  different = lo-same
  multiple = dim(peaks_lo_multiple)[1] #appeared in liftover results on multiple chromosomes
  failed = total - sum(lo,multiple) #not showing up in liftover results at all
  Category <- c('Same width','Different widths','Multiple','Failed','Total')
  Count <- c(same,different,multiple,failed,total)
  summary <- data.frame(Category,Count)
  print(paste0('Peak half-width: ', half_width))
  if (show) {print(summary)}
  #Return granges objects for each category
  names_same = rownames(df_same)
  names_dif = rownames(df_dif)
  names_multiple = peaks_lo_multiple$name
  names_failed_multiple = setdiff(peaks$name,peaks_lo$name)
  names_failed = setdiff(names_failed_multiple,names_multiple)
  return(list(summary = summary, names_same = names_same, names_dif = names_dif, names_multiple = names_multiple, names_failed=names_failed))
}

#Summarize where peaks were lost in round 1, part 2 - lifting back to individual species genomes
summarizeMultiConsensusPeaksPart2 <- function(species,consensus_peaks,consensus_df, species_con_summits1,species_con_summits2,
                                              species_con_peaks1,species_con_peaks2, species_con_peaks){
  df_part2_sp = list()
  exclusion_list = list()
  df_total_failed= list()
  for (sp in species){
    df_stepB = length(consensus_peaks) - length(species_con_summits1[[sp]]) #peaks that failed liftover (failed and multiple)
    df_step4a = length(species_con_summits1[[sp]]) - length(species_con_summits2[[sp]]) # peaks that failed liftover (different widths)
    df_step4b = length(species_con_summits2[[sp]]) - length(species_con_peaks1[[sp]]) # peaks that failed liftover (overlapping ends of chromosomes)
    problem_indices = which(consensus_df$liftback_result=='Problem'&consensus_df$liftback_species==sp)
    problem_names = consensus_df$name[problem_indices]
    df_step5_problem = length(unique(problem_names))
    blacklist_indices = which(consensus_df$liftback_result=='Blacklist'&consensus_df$liftback_species==sp)
    blacklist_names = consensus_df$name[blacklist_indices]
    df_step5_blacklist = length(unique(blacklist_names))
    chr_indices = which(consensus_df$liftback_result=='Chr_removed'&consensus_df$liftback_species==sp)
    chr_names = consensus_df$name[chr_indices]
    df_step5_chr = length(unique(chr_names))
    unique_step5 = length(unique(c(problem_names, blacklist_names, chr_names)))
    exclusion_indices = which(consensus_df$liftback_result !='' & consensus_df$liftback_species==sp)
    exclusion_names = consensus_df$name[exclusion_indices]
    exclusion_list[[sp]] = unique(exclusion_names)
    df_total_failed[[sp]] = length(consensus_peaks)- length(species_con_peaks[[sp]])
    df_part2_sp[[sp]] = c(length(consensus_peaks),df_stepB,df_step4a,df_step4b,df_step5_problem, 
                          df_step5_blacklist, df_step5_chr, length(exclusion_list[[sp]]),
                          df_total_failed[[sp]],length(species_con_peaks[[sp]]),
                          length(species_con_peaks[[sp]])/length(consensus_peaks))
    #check that exclusion lists for each species overlap by correct amount so that removed for each species - length of exclusion list = total removed
    total_all_steps = df_stepB +df_step4a + df_step4b + unique_step5
    if (total_all_steps==length(exclusion_list[[sp]])) {
      print(paste0('Part 2 1st check successful - total peaks removed in species ', sp, ' is equal to the sum of peaks removed from each step.'))
    } else {warning(paste0('Part 2 1st check failed - total peaks removed in species ', sp, ' is NOT equal to the sum of peaks removed from each step for both species.'))}
  }
  df_part2 <- do.call(cbind, df_part2_sp)
  rownames(df_part2) = c('Total to start','Removed, Step B','Removed, Step 4a','Removed, Step 4b', 'Removed, Step 5 (problem)',
                         'Removed, Step 5 (blacklist)','Removed, Step 5 (non-standard chromosomes)', 
                         'Total removed (each species)','Total removed (all species)','Total remaining','Proportion remaining')
  unique_exclusion = unique(unlist(exclusion_list))
  all_same_length <- length(unique(unlist(df_total_failed)))==1
  if (df_total_failed[[1]]== length(unique_exclusion) & all_same_length){
    print(paste0('Part 2 2nd check successful - total peaks removed for all species is equal to uniquely excluded peaks.'))
  } else {warning(paste0('Part 2 2nd check failed - total peaks removed per species is NOT equal to uniquely excluded peaks.'))}
  return(df_part2)
}


#Updates metadata to show which consensus peaks originated from which species
updateMetadataPostMerge <- function(mergedPeaks, metadataList, originalGrangesList, speciesNames) {
  cleanedMetadataList <- list()
  
  # Iterate over species names
  for (speciesIdx in seq_along(speciesNames)) {
    species = speciesNames[speciesIdx]
    metadata = metadataList[[species]]
    
    # Extract all original names from mergedPeaks and match them with the species metadata
    splitOriginalNames <- lapply(mcols(mergedPeaks)$orig_name, function(name) strsplit(name, ";")[[1]])
    originalNames <- unlist(splitOriginalNames)  # Flatten the list of original names
    indices <- match(originalNames, metadata$original_peak_name)  # Find indices in metadata
    valid_indices <- !is.na(indices)
    
    # Create a dataframe to collect matched metadata entries
    if (any(valid_indices)) {
      relevant_metadata = metadata[indices[valid_indices], , drop = FALSE]
      valid_peak_indices <- rep(seq_along(splitOriginalNames), sapply(splitOriginalNames, length))[valid_indices]
      
      # Calculate original summits vectorized
      originalSummits = (relevant_metadata$original_start + relevant_metadata$original_end) / 2
      
      # Calculate new summits vectorized
      mergedStarts = start(mergedPeaks[valid_peak_indices])
      mergedEnds = end(mergedPeaks[valid_peak_indices])
      newSummits = (mergedStarts + mergedEnds) / 2
      
      # Calculate summit adjustments vectorized
      summitAdjustments = newSummits - originalSummits
      
      # Apply strand flip adjustments in a vectorized way
      strandFlips = relevant_metadata$strand_flip
      summitAdjustments[strandFlips] = -summitAdjustments[strandFlips]
      
      # Update relevant metadata
      relevant_metadata$summit_adj = summitAdjustments
      relevant_metadata$consensus_name = names(mergedPeaks)[valid_peak_indices]
      
      # Store cleaned metadata in the list
      cleanedMetadataList[[species]] = relevant_metadata
    } else {
      cleanedMetadataList[[species]] = metadata[FALSE, ]  # Empty metadata frame with the same structure
    }
    print(paste0('Updated metadata for ',speciesNames[[speciesIdx]]))
  }
  
  return(cleanedMetadataList)
}