#Folders----
folder = 'Example_dataset/' #should be the same as current_folder in run_CrossPeak_example.sh
figure_folder = paste0(folder,'Figures/'); dir.create(figure_folder)
output_folder = paste0(folder,'Output/'); dir.create(output_folder) #should be the same as output_folder in run_CrossPeak_example.sh

#Parameters----
# *May need to be modified, depending on your analysis* See README for full descriptions.
peaks_format = 'bed' #'bed' or 'rds' - filename extension of original peaks
summit_half_width = 5 #Half-width for summits so 11-bp summits will be used for round 1 liftovers
peak_half_width = 250 #Half-width for peaks so peaks are 501-bp with summits exactly centered
indel_tol = 5 #5-bp tolerance for summit adjustment during liftovers
overlap_tol = 300 #Set to 0 if you do not want any overlapping peaks in the final set
remove_non_standard_chr = TRUE #Will remove any chromosomes that do not start with "chr" or that start with "chr" but then also include an underscore
chr_to_remove = c('chrM', 'chrY') #Will remove all peaks from any species on these chromosomes
summit_half_width_rd2 = 25 #51-bp summits are used for round 2 liftovers 
indel_low = 5; indel_med = 50; indel_high = 450 #Indel tolerances for round 2 (peaks meeting any of these tolerances will be added to the consensus set and the category will be saved in metadata)

#Setup species----
# *May need to be modified, depending on your analysis* See README for full details
species = c('human','chimp','rhesus') #For more than 2 species, should be in order of evolutionary distance (but order of first 2 does not matter)
prefixes = c('Hu', 'Ch', 'Rh'); names(prefixes) <- species 
genome_names <- c(human = 'Human', chimp = 'Chimp', rhesus = 'Rhesus')#names of genomes for granges objects (can be official genome build names or names from Haltoolkit)
genome_ancestral = 'primate' #name of ancestral genome in hal file
