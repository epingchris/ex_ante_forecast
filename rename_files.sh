#!/bin/bash

# Define the base directory
BASE_DIR="/maps/epr26/tmf_pipe_out"

# Loop through each subdirectory in the base directory
for dir in "$BASE_DIR"/*/; do
  # Extract the folder name (e.g., '1047')
  folder_name=$(basename "$dir")

  # Loop through all files in the current subdirectory
  for file in "$dir"/*; do
    # Get the base name of the file (e.g., '1047additionality.csv')
    base_file_name=$(basename "$file")
    
    # Check if the file starts with the folder name
    if [[ "$base_file_name" == "$folder_name"* ]]; then
      # Generate the new file name by removing the folder name prefix
      new_file_name="${base_file_name#$folder_name}"
      
      # Rename the file
      mv "$file" "$dir$new_file_name"
      echo "Renamed '$file' to '$dir$new_file_name'"
    fi
  done
done
