#!/bin/bash

# Setup species, folders, and files *May need to modify for your analysis*
script_folder=$(pwd)         # Base folder where the R scripts are located
folder="./Example_dataset" # Folder to save all intermediate results and logs
output_folder="./Example_dataset/Output" # Folder to save results
params_file="CrossPeak_parameters.R" # Leave unchanged to use the default parameters file, or you can create your own parameters file
HAL_FILE="${folder}/primates.hal" # HAL file with species' genomes and ancestral genome
species=("human" "chimp" "rhesus") # Lowercase versions of names of species in HAL file
ANCESTRAL_GENOME="primate" # Name of ancestral genome in HAL file
docker_image="quay.io/comparative-genomics-toolkit/cactus:v2.9.3" # Docker image for HAL tools
use_docker=true # Set to true to use Docker, false to use local HAL installation
cleanup=true #Set to true to clean up intermediate files after the pipeline finishes. Set to false to assist in debugging

# Check for Docker if using Docker
if [[ "$use_docker" == true ]]; then
  if ! command -v docker &> /dev/null; then
    echo "Error: Docker is not installed or not in the PATH. Please install Docker or set use_docker=false."
    exit 1
  fi
fi

# Check for HAL tools if not using Docker
if [[ "$use_docker" == false ]]; then
  if ! command -v halLiftover &> /dev/null; then
    echo "Error: halLiftover is not installed or not in the PATH. Please install HAL tools or set use_docker=true."
    exit 1
  fi
fi

# Create the output folder if it doesn't exist
if [[ ! -d "$folder" ]]; then
  mkdir -p "$folder"
  echo "Created output folder: $folder"
fi

# Record the files in the current folder at the beginning
initial_files=$(mktemp) # Temporary file to save the list of initial files
ls "$folder" > "$initial_files"

# Step 1: Run R script to generate initial data
echo "Running Step 1 in R..."
Rscript "${script_folder}/CrossPeak_Step1.R" "$script_folder" "$params_file" > "${folder}/step1.log" 2>&1
if [[ $? -ne 0 ]]; then
  echo "Error: Step 1 failed. Check ${folder}/step1.log for details."
  exit 1
fi

# Step A: Run halLiftover for each species
echo "Running Step A: halLiftover..."
hal_liftover_log="${folder}/stepA_halLiftover.log"
: > "$hal_liftover_log"
for sp in "${species[@]}"; do
  input_file="${folder}/${sp}_summits.bed"
  output_file="${folder}/${sp}_summits_lo.bed"

  # Check if input file exists before running halLiftover
  if [[ ! -f "$input_file" ]]; then
    echo "Error: Input file $input_file not found for species $sp."
    exit 1
  fi

  # Run halLiftover, either with Docker or directly
  echo "Lifting over $sp..." | tee -a "$hal_liftover_log"
  if [[ "$use_docker" == true ]]; then
    docker run --rm -v "$folder:$folder" -w "$script_folder" "$docker_image" halLiftover "$HAL_FILE" "${sp^}" "$input_file" "$ANCESTRAL_GENOME" "$output_file"  >> "$hal_liftover_log" 2>&1
  else
    halLiftover "$HAL_FILE" "${sp^}" "$input_file" "$ANCESTRAL_GENOME" "$output_file" >> "$hal_liftover_log" 2>&1
  fi

  # Check if halLiftover was successful or the output BED file is empty
  if [[ $? -ne 0 || ! -s "$output_file" ]]; then
    echo "Error: halLiftover failed for $sp, or the output file $output_file is empty. Check inputs and HAL file."
    exit 1
  fi
done

# Steps 2 and 3: Run the next R script
echo "Running Steps 2 and 3 in R..."
Rscript "${script_folder}/CrossPeak_Step2and3.R" "$script_folder" "$params_file" > "${folder}/step2and3.log" 2>&1
if [[ $? -ne 0 ]]; then
  echo "Error: Steps 2 and 3 failed. Check ${folder}/step2and3.log for details."
  exit 1
fi

# Step B: Run halLiftover for consensus summits
echo "Running Step B: halLiftover for consensus summits..."
hal_liftover_log="${folder}/stepB_halLiftover.log"
: > "$hal_liftover_log"
consensus_input_file="${folder}/consensus_summits.bed"
if [[ ! -f "$consensus_input_file" ]]; then
  echo "Error: Input file $consensus_input_file not found. Ensure Steps 2 and 3 completed successfully."
  exit 1
fi

for sp in "${species[@]}"; do
  consensus_output_file="${folder}/${sp}_consensus_summits_lo.bed"

  # Run halLiftover, either with Docker or directly
  echo "Lifting over consensus summits for $sp..." | tee -a "$hal_liftover_log"
  if [[ "$use_docker" == true ]]; then
    docker run --rm -v "$folder:$folder" -w "$script_folder" "$docker_image" halLiftover "$HAL_FILE" "$ANCESTRAL_GENOME" "$consensus_input_file" "${sp^}" "$consensus_output_file"  >> "$hal_liftover_log" 2>&1
  else
    halLiftover "$HAL_FILE" "$ANCESTRAL_GENOME" "$consensus_input_file" "${sp^}" "$consensus_output_file"  >> "$hal_liftover_log" 2>&1
  fi

  # Check if halLiftover was successful or the output BED file is empty
  if [[ $? -ne 0 || ! -s "$consensus_output_file" ]]; then
    echo "Error: halLiftover failed for consensus summits for $sp, or the output file $consensus_output_file is empty. Check inputs and HAL file."
    exit 1
  fi
done

# Steps 4, 5 and 6: Run the next R script
echo "Running Steps 4, 5, and 6 in R..."
Rscript "${script_folder}/CrossPeak_Step4and5and6.R" "$script_folder" "$params_file" > "${folder}/step4and5and6.log" 2>&1
if [[ $? -ne 0 ]]; then
  echo "Error: Steps 4, 5, or 6 failed. Check ${folder}/step4and5and6.log for details."
  exit 1
fi

# Step C: Run halLiftover for species-specific peaks and summits
echo "Running Step C: halLiftover for species-specific peaks and summits..."
hal_liftover_log="${folder}/stepC_halLiftover.log"
: > "$hal_liftover_log"
for source_species in "${species[@]}"; do
  for target_species in "${species[@]}"; do
    # Skip if source and target species are the same
    if [[ "$source_species" == "$target_species" ]]; then
      continue
    fi

    # Generate file names for summits
    input_summits_file="${folder}/${source_species}_spec_summits.bed"
    output_summits_file="${folder}/${source_species}_on_${target_species}_spec_summits_lo.bed"

    # Check if input file exists for summits
    if [[ ! -f "$input_summits_file" ]]; then
      echo "Error: Input file $input_summits_file not found for $source_species." | tee -a "$hal_liftover_log"
      exit 1
    fi

    # Run halLiftover for summits
    echo "Lifting over summits from $source_species to $target_species..." | tee -a "$hal_liftover_log"
    if [[ "$use_docker" == true ]]; then
      docker run --rm -v "$folder:$folder" -w "$script_folder" "$docker_image" halLiftover "$HAL_FILE" "${source_species^}" "$input_summits_file" "${target_species^}" "$output_summits_file" >> "$hal_liftover_log" 2>&1
    else
      halLiftover "$HAL_FILE" "${source_species^}" "$input_summits_file" "${target_species^}" "$output_summits_file" >> "$hal_liftover_log" 2>&1
    fi

    # Check if halLiftover was successful or the output BED file is empty
    if [[ $? -ne 0 || ! -s "$output_summits_file" ]]; then
      echo "Error: halLiftover failed or output file $output_summits_file is empty. Check inputs and HAL file." | tee -a "$hal_liftover_log"
      exit 1
    fi

    # Generate file names for peaks
    input_peaks_file="${folder}/${source_species}_spec_peaks.bed"
    output_peaks_file="${folder}/${source_species}_on_${target_species}_spec_peaks_lo.bed"

    # Check if input file exists for peaks
    if [[ ! -f "$input_peaks_file" ]]; then
      echo "Error: Input file $input_peaks_file not found for $source_species." | tee -a "$hal_liftover_log"
      exit 1
    fi

    # Run halLiftover for peaks
    echo "Lifting over peaks from $source_species to $target_species..." | tee -a "$hal_liftover_log"
    if [[ "$use_docker" == true ]]; then
      docker run --rm -v "$folder:$folder" -w "$script_folder" "$docker_image" halLiftover "$HAL_FILE" "${source_species^}" "$input_peaks_file" "${target_species^}" "$output_peaks_file" >> "$hal_liftover_log" 2>&1
    else
      halLiftover "$HAL_FILE" "${source_species^}" "$input_peaks_file" "${target_species^}" "$output_peaks_file" >> "$hal_liftover_log" 2>&1
    fi

    # Check if halLiftover was successful or the output BED file is empty
    if [[ $? -ne 0 || ! -s "$output_peaks_file" ]]; then
      echo "Error: halLiftover failed or output file $output_peaks_file is empty. Check inputs and HAL file." | tee -a "$hal_liftover_log"
      exit 1
    fi
  done
done

# Step 7: Run the next R script
echo "Running Step 7 in R..."
Rscript "${script_folder}/CrossPeak_Step7.R" "$script_folder" "$params_file" > "${folder}/step7.log" 2>&1
if [[ $? -ne 0 ]]; then
  echo "Error: Step 7 failed. Check ${folder}/step7.log for details."
  exit 1
fi

# Step D: Reciprocal liftover for round 2 peaks
echo "Running Step D: Reciprocal halLiftover for round 2 peaks..."
hal_liftover_log="${folder}/stepD_halLiftover.log"
: > "$hal_liftover_log"
for source_species in "${species[@]}"; do
  for target_species in "${species[@]}"; do
    # Skip if source and target species are the same
    if [[ "$source_species" == "$target_species" ]]; then
      continue
    fi

    # Generate file names
    input_file="${folder}/${target_species}_on_${source_species}_rd2_con_peaks.bed"
    output_file="${folder}/${target_species}_from_${source_species}_rd2_con_peaks_lo.bed"

    # Check if input file exists
    if [[ ! -f "$input_file" ]]; then
      echo "Error: Input file $input_file not found for $source_species to $target_species." | tee -a "$hal_liftover_log"
      exit 1
    fi

    # Run halLiftover, either with Docker or directly
    echo "Lifting over round 2 consensus peaks from $target_species to $source_species..." | tee -a "$hal_liftover_log"
    if [[ "$use_docker" == true ]]; then
      docker run --rm -v "$folder:$folder" -w "$script_folder" "$docker_image" halLiftover "$HAL_FILE" "${source_species^}" "$input_file" "${target_species^}" "$output_file" >> "$hal_liftover_log" 2>&1
    else
      halLiftover "$HAL_FILE" "${source_species^}" "$input_file" "${target_species^}" "$output_file" >> "$hal_liftover_log" 2>&1
    fi

    # Check if halLiftover was successful or the output BED file is empty
    if [[ $? -ne 0 || ! -s "$output_file" ]]; then
      echo "Error: halLiftover failed for $source_species to $target_species, or the output file $output_file is empty. Check inputs and HAL file." | tee -a "$hal_liftover_log"
      exit 1
    fi
  done
done

# Steps 8 and 9: Run the next R script
echo "Running Steps 8 and 9 in R..."
Rscript "${script_folder}/CrossPeak_Step8and9.R" "$script_folder" "$params_file" > "${folder}/step8and9.log" 2>&1
if [[ $? -ne 0 ]]; then
  echo "Error: Step 8 or 9 failed. Check ${folder}/step8and9.log for details."
  exit 1
fi


# Cleanup: Delete intermediate files if the cleanup flag is set
if [[ "$cleanup" == true ]]; then
  echo "Performing cleanup of intermediate files..."

  # Get the list of .log files to keep
  log_files=$(find "$folder" -maxdepth 1 -type f -name "*.log")

  # Get the list of files to keep (initial files + .log files)
  files_to_keep=$(mktemp)
  cat "$initial_files" > "$files_to_keep"
  echo "$log_files" >> "$files_to_keep"

  # Find and remove all files in the current folder but not output folder except those in the "files_to_keep"
  for file in $(find "$folder" -maxdepth 1 -type f); do
    if ! grep -qF "$(basename "$file")" "$files_to_keep"; then
      rm -f "$file"
    fi
  done

  # Clean up temporary files
  rm -f "$initial_files" "$files_to_keep"

  echo "Cleanup completed. All intermediate files removed."
else
  echo "Cleanup flag is set to false. Intermediate files are retained in ${folder}."
fi

# Final message
echo "CrossPeak completed successfully. Logs are in ${folder} and output files are in ${output_folder}."
