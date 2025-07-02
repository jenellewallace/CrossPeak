# CrossPeak

## Overview

CrossPeak is a tool to generate cross-species consensus ATAC-seq peak sets for 2 or more species.

## Installation

Start R and clone the CrossPeak repository:
```R
devtools::install_github("jenellewallace/CrossPeak@main")
````
CrossPeak uses renv for package management and reproducibility. Install renv (https://github.com/rstudio/renv) if you have not already and use R to install all dependencies:

```R
install.packages('renv')
renv::restore()
````
CrossPeak requires the creation of an inferred ancestral genome and a way to lift over coordinates between individual species genomes and the ancestral genome. Although not strictly required by CrossPeak, our recommended workflow uses Hierarchical Alignment (HAL) format to create the ancestral genome and halLiftover to lift over coordinates. If you would like to use HAL with CrossPeak, follow the HAL installation instructions at https://github.com/ComparativeGenomicsToolkit/hal. If you run into any problems with installing HAL, we recommend you instead install HAL as part of cactus (https://github.com/ComparativeGenomicsToolkit/cactus/releases) via a Docker image. We provide an example script that runs both the R steps and the halLiftover steps in one script and includes options for running the halLiftover steps with a local HAL installation or via the Docker image.

## Citing

CrossPeak is described in the following preprint:
Interspecies Organoids Reveal Human-Specific Molecular Features of Dopaminergic Neuron Development and Vulnerability
Sara Nolbrant, Jenelle L. Wallace, Jingwen Ding, Tianjia Zhu, Jess L. Sevetson, Janko Kajtez, Isabella A. Baldacci, Emily K. Corrigan, Kaylynn Hoglin, Reed McMullen, Matthew T. Schmitz, Arnar Breevoort, Dani Swope, Fengxia Wu, Bryan J. Pavlovic, Sofie R. Salama, Agnete Kirkeby, Hao Huang, Nathan K. Schaefer, Alex A. Pollen
bioRxiv 2024.11.14.623592; doi: https://doi.org/10.1101/2024.11.14.623592

If you use HAL tools (https://github.com/ComparativeGenomicsToolkit/hal) with CrossPeak, then also cite:
Glenn Hickey, Benedict Paten, Dent Earl, Daniel Zerbino, and David Haussler. HAL: A Hierarchical Format for Storing and Analyzing Multiple Genome Alignments. Bioinformatics. 2013. https://doi.org/10.1093/bioinformatics/btt128

If you use the primate HAL file we provide in the Example_dataset folder for your own analysis, then also cite the papers describing those genome alignments:

Mao, Y., Catacchio, C.R., Hillier, L.W. et al. A high-quality bonobo genome refines the analysis of hominid evolution. Nature 594, 77–81 (2021). https://doi.org/10.1038/s41586-021-03519-x

Zev N. Kronenberg et al. ,High-resolution comparative analysis of great ape genomes. Science360,eaar6343(2018). https://doi.org/10.1126/science.aar6343

## Running CrossPeak

### Running the example script:
We provide an example Shell script that runs all of the CrossPeak and halLiftover steps with a single command: run_CrossPeak_example.sh. To run this script on the example dataset, users will need install CrossPeak and HAL and download the primate HAL file into the Example_datasets folder: https://www.dropbox.com/scl/fi/xp97bzed9tilmaluvllog/primates.hal?rlkey=vknl7zg0t36m3eb2d8d5tfo5c&st=z8gnao1x&dl=0 (note this file is ~30GB and may take a while to download). The file CrossPeak_parameters.R includes default parameters that we use to run the pipeline in our own hands so we recommend leaving this file unchanged for running the example dataset for the first time.

### Running CrossPeak on your own dataset:
 To run CrossPeak on your own dataset, you will need to modify the files run_CrossPeak_example.sh and CrossPeak_parameters.R to include the species in your analysis and your chosen parameters. The steps described below correspond to the steps listed in example script. Numbered steps take place in R while lettered steps (liftovers) require steps on the command line using HAL tools or other software. Steps 1-5 represent round 1 and steps 6-7 represent round 2 (described below). Numbered steps after a lettered step (for example, Step 2 which occurs after Step A) require bed files as input (specifics listed for each step below). If you use our recommended workflow with HAL and our example Shell script, correctly named files should be seamlessly imported and exported from R. If you use a different method for liftovers, then you will need to run the numbered R scripts individually and perform the lfitovers for lettered steps in the order described below, naming files accordingly.

## Method details
Summary of the CrossPeak pipeline:
![CrossPeak Summary](README_figure/CrossPeak_summary.png)

### Preprocessing:
CrossPeak takes as input a single peak set of summit-centered ATAC-seq peaks for each species in the analysis. The peak set can be generated with any tool, but the summits must be exactly centered with equal numbers of base pairs on either side (for example, peaks could be 501 bp in width with a 250-bp half-width). If using MACS2 with standard parameters, the peaks will need to be postprocessed to center the summits and chose a fixed width. If using CrossPeak with single cell ATAC-seq data with peaks called for each cell type, the user must choose whether they want one set of cross-species consensus peaks for each cell type (in which case the pipeline must be run once for each cell type) or whether they want one set of cross-species consensus peaks for the entire dataset, in which case they must collapse the peak sets for each cell type into one unified peak set across all cell types. We suggest using iterative overlap peak merging (https://www.archrproject.com/bookdown/the-iterative-overlap-peak-merging-procedure.html) for this purpose, which is implemented in packages like ArchR and ScenicPlus.

#### Setup: Choose species and save data structures to working folder
CrossPeak can work with any species for which you have an inferred ancestral genome, a way to lift over coordinates between individual species genomes and the ancestral genome, and a blacklist of regions to exclude. We currently provide bed files with blacklists for human, chimpanzee, and macaque. To modify the pipeline to work with other species:
- Update section called “Setup species and blacklists” in this script
- If you would like to exclude peaks on non-standard chromosomes, update the function called findNonStandardChromosomes if applicable (current version works well for primate genomes)
Species should be provided in order of evolutionary relatedness with the species of greatest interest first and outgroups last because CrossPeak uses the species order to prioritize peak merging (see Step 3 below). For example, if we are most interested in identifying human-specific peaks, we would use the order: human, chimpanzee, macaque. 

#### Step 1: Import peaks files and make summits
Import: [species]_peaks.rds or [species]_peaks.bed, Export: [species]_summits.bed
Peak files may be in .bed format or GenomicRanges objects in .rds format. Peak files are imported, peaks overlapping blacklisted regions or located on user-excluded chromosomes are excluded, and a summits file is created and exported. Summit regions consist of a shorter region around the single base pair summit that will be preserved across species. In preliminary tests, we found that shorter summit region widths offered the best compromise between precise summit localization and liftover precision (success of reciprocal liftovers and minimization of liftovers that failed or went to multiple locations), so we recommend using 11-bp summits (this parameter can be adjusted by the user).

#### Step A: Liftover from individual species genomes to ancestral genome
Summits are lifted over from each individual species genome to the same ancestral genome.

#### Step 2: Concatenate regions and make summits from peaks
Import: [species]_summits_lo.bed
We take the union of the lifted ranges and exclude peaks that failed liftover, lifted to multiple chromosomes, or whose lifted summit widths are more than an indel tolerance (indel_tol, default is 5 bp) different from the original summit widths. Then, we rebuild peaks on the ancestral genome by extending the summits by 250 bp in each direction.

#### Step 3: Merge peaks on ancestral genome to make the consensus set
Export: consensus_summits_lo.bed
To explore parameters for peak merging, we suggest plotting histograms of the peak overlap widths for each pair of species and choosing a value for the overlap tolerance that marks the transition from a nearly uniform to a nonlinearly increasing distribution. We reason that peaks with less overlap than this value likely represent nearby but distinct peaks in each species whereas peaks with more overlap are more likely to be the same peak in each species with the summits slightly offset or mislocalized. The optimal value for the overlap tolerance (overlap_tol) may vary depending on the dataset and scientific goals. We have chosen values between 300 and 350 bp for our datasets, meaning that peaks with greater overlap than this value will be merged and peaks with less overlap will be retained. This value can be set to 0 if the user prefers not to have any overlapping peaks in the final set. 
CrossPeak merges peaks in order of evolutionary relatedness, given by the order of species that was input by the user. In our example, we first merge human and chimpanzee peaks that overlap by more than the overlap tolerance by taking the average summit position. Macaque peaks are added to the consensus set only if they do not overlap human-chimpanzee peaks by more than the overlap tolerance and so on until all species peaks have been processed. After peak merging is complete, consensus summit regions are created and exported.

#### Step B: Liftover the consensus peaks from the ancestral genome to individual species genomes
Consensus summit regions are lifted over from the ancestral genome to each individual species genome.

#### Step 4: Concatenate regions and make summits from peaks
Import: [species]_consensus_summits_[overlap_tol]_lo.bed
We take the union of the lifted ranges and exclude peaks that failed liftover, lifted to multiple chromosomes, or whose lifted summit widths are more than an indel tolerance (default is 5 bp) different from the original summit widths. Then, we rebuild peaks on the individual species genomes by extending the summits by 250 bp in each direction.

#### Step 5: Remove peaks that do not meet criteria from both species
We compare the location of the consensus summits with the original summit location and exclude peaks with summits that are too far from the original summit or overlap blacklisted regions in any species. After this step, round 1 is complete and peaks in the consensus set are “round 1 consensus peaks.” See description under Outputs below. If you only care about the highest confidence consensus peaks and do not care about species-specific peaks, you could end the analysis here.

#### Step 6: Combine potentially species-specific peaks that failed to make it into consensus set
Export: [species]_spec_summits.bed
For each species, we combine peaks that failed to lift from that species’ genome to the ancestral genome in Step A and peaks that originated from that species’ genome, lifted to the ancestral genome, but failed to lift back to one of the other species’ genomes in Step B. These represent potentially species-specific peaks that will be explored further.

#### Step C: Liftover from each species genome directly to other species’ genomes
Import: [species]_on_[other species]_spec_summits_lo.bed, [species]_on_[other species]_spec_peaks_lo.bed
Longer summits (we recommend 51-bp) as well as the full-size peaks are lifted over directly between each pair of species genomes. 

#### Step 7: Divide peaks from step 6 into ones that can be added to consensus set (with variable tolerances tagged separately) and species-only peaks
In order to include as many peaks as possible in the consensus set, we offer the option to relax the indel tolerance in this step and to retain peaks that lift over to other species’ genomes but with less precise summit localization than round 1 peaks. We tag these peaks with separate metadata categories. Available categories with parameters that can be adjusted by the user are:
1. round2_indel_low: We recommend setting indel_low to the same indel_tol as used for round 1. This category represents peaks that lifted within the same indel tolerance as round 1 peaks (suggests an issue with the ancestral genome in this location but the peak location is conserved).
2. round2_indel_med: Peaks that lifted with an intermediate indel tolerance (set to 50 bp in our example).
3. round2_indel_high: Peaks lifted with with a high indel tolerance (set to 450 bp in our example).
4. round2_summit_deletion:  Peaks where the summit failed but the whole peak was lifted over successfully
5. round2_failed: Peaks that failed to lift over entirely.
6. round2_summit_with_large_indel: Peaks where the summit region lifted but the resulting width was so large that the summit would move outside of the existing peak width  
7. round2_summit_failed_large_indel_peaks: Peaks where the summit region failed to lift and the width of the lifted peak was so large that the summit would move outside of the existing peak width
Peaks that fall in categories 1-4 for all pairs of species are added to the consensus peak set with the category for each species included in the metadata. Peaks that fall in categories 5-7 for all pairs of species are classified as “species-only peaks” for that species. Peaks that fall into different sets of categories (for example, category 1 for one species and category 6 for another species) are excluded as are peaks with multi-mapping summits or failed summits and multi-mapping peaks. Round 2 peaks that overlap round 1 peaks by more than the overlap tolerance are also excluded. Finally, peaks that overlap blacklisted regions or are located on excluded chromosomes for any species are excluded for all species.

#### Step D: Reciprocal liftover for round 2 peaks
Lifts the results of step 7 (round 2 consensus peaks) back to the original species genome (reciprocal liftover confirms a 1:1 relationship across species for all consensus peaks).


#### Step 8: Final cleanup and export results
Similar to Step 5, this step removes any round 2 peaks that failed to lift back to the correct location on the original species' genome. It also removes any round 2 peaks that overlap round 1 peaks by more than the overlap tolerance.

#### Step 9: Prints summaries of how many peaks were retained and how many peaks were lost at each step

#### Outputs:
Files are saved in the Output folder as GenomicRanges objects in .rds files with the following names:
```R
species_rd1_peaks.rds
species_allcon_peaks.rds
species_only_peaks.rds
````
Bed files for each species are also exported to the same folder.
Intermediate files are saved at each step (with the option to clean up or leave these files at the end) but the final outputs are:
species_rd1_peaks: Consensus peaks from round 1 only in the form of a list with slots for each species that contains all round 1 and round 2 consensus peaks in the coordinates of each individual species genomes with corresponding names (i.e. Peak_1 has the same name and is the corresponding peak for all species). Therefore, all slots will have the exact same number of peaks.
species_allcon_peaks: Same as above but including all round 1 and round 2 consensus peaks.
Metadata columns are:
name: consensus name for the peak that corresponds across species (Peak_1, etc). Note that peak numbers are unique but will not always be consecutive since peaks are numbered in step 3 but some will be excluded in steps 4-5.
orig_name: original name corresponding to the names in the peak files for each species created in step 1 (e.g. Hu_peak_1, Ch_peak_1, etc). More than one name may be listed here, separated by semicolons, if the consensus peak involved merging more than one original peak.
category: category referring to the lift over status – will be “round1_indel_low” if the peak was from round 1 and otherwise it will be the round 2 categories assigned in step 7 (one category for the liftover from the given species to each other species in the original order provided, separated by semicolons).
Species_only_peaks: Species-only peaks in the form of a list with slots for each species that contain species-only peaks for each species. These are the peaks that could not be lifted over successfully to the other species’ genomes; therefore, each slot may contain a different number of peaks for each species. Metadata columns are:
name: final name for the peak. The names here will be unique compared to the consensus names and compared to the other slots. For example, if the consensus peaks are Peak_1 to Peak_100, then the human-only peaks might be Peak_101 to Peak_200, and the chimpanzee-only peaks would be Peak_201 to Peak_300, etc.
orig_name: original name corresponding to the names in the peak files for each species created in step 1. All the names in this slot should indicate that the peak originated from that species (for example species_only_peaks$human should only include original names that begin with Hu_peak).
category: category referring to the lift over status. All species-only peaks should have categories #5,6,or 7 from the list in Step 7. 

### Downstream analysis:
For differential accessibility analysis, after running CrossPeak, users will need to quantify counts for the consensus peaks and species-only peak sets for each cell/sample to obtain a counts matrix. For example, for single cell ATAC-seq data, this can be performed with Signac FeatureMatrix or other tools. Then users can perform differential accessibility analysis on the consensus peak set using existing tools to test which peaks have increased or decreased accessibility across species, among other analyses. See our preprint listed above for examples of downstream analysis using CrossPeak outputs.

