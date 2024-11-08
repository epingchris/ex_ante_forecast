import os
import re
import pandas as pd
import numpy as np
from glob import glob
from datetime import datetime
import geojson
from area import area

# Function to find new projects to process (those not already in the CSV)
def find_new_projects(directories):
    # Get list of monitored project names in the DataFrame
    existing_projects = set(df.loc[df['build_m_table_memory_GB'].notnull(), 'proj_name'].tolist())
    new_project_paths = set()
    
    # Get list of finished project names based on subfolders in directory, excluding subfolders that include old results or special data
    project_subfolders = []
    project_directories = []
    for directory in directories:
        for entry in os.scandir(directory):
            if entry.is_dir():
                name = entry.name
                if re.match(r"as\d{2}$", name):  # Exclude "as[0-9][0-9]"
                    continue
                if name.startswith(("slopes", "rescaled", "srtm")):  # Exclude "slopes" and "rescaled" and "srtm"
                    continue
                project_subfolders.append(name)
                project_directories.append(directory)
    
    for project_subfolder, project_directory in zip(project_subfolders, project_directories):
        if not project_subfolder in existing_projects:
            new_project_paths.add(os.path.join(project_directory, project_subfolder))
            
    return new_project_paths

# Function to parse memory usage and runtime from log files
def parse_log_file(log_file_path):
    memory_usage = {
        'calculate_k_memory_GB': np.nan,
        'find_potential_matches_memory_GB': np.nan,
        'build_m_table_memory_GB': np.nan,
        'find_pairs_memory_GB': np.nan
    }
    execution_time_seconds = np.nan
    
    with open(log_file_path, 'r') as file:
        for line in file:
            # Extract memory usage information
            if line.startswith("Memory used by calculate_k"):
                memory_usage['calculate_k_memory_GB'] = float(line.split()[-2]) / 1e6  # Convert KB to GB
            elif line.startswith("Memory used by find_potential_matches"):
                memory_usage['find_potential_matches_memory_GB'] = float(line.split()[-2]) / 1e6
            elif line.startswith("Memory used by build_m_table"):
                memory_usage['build_m_table_memory_GB'] = float(line.split()[-2]) / 1e6
            elif line.startswith("Memory used by find_pairs"):
                memory_usage['find_pairs_memory_GB'] = float(line.split()[-2]) / 1e6
            
            # Extract execution time information
            elif line.strip().endswith("seconds."):
                if "Total execution time:" in line:
                    execution_time_seconds = float(line.split()[-2])
    
    return memory_usage, execution_time_seconds

# Function to find the latest log file based on date
def find_latest_log_file(project_name, folder_path):
    # Define the pattern for the log files we want to find
    pattern = f"out_{project_name}_*_out.txt"
    log_files = glob(os.path.join(folder_path, pattern))
    
    # Check if any files match the pattern
    if not log_files:
        print("in find_latest_log_file function: log_files is empty")
        return None
    
    # Extract date from each file name and find the latest one
    date_format = "%Y_%m_%d"  # Assuming date in filename is in format YYYYMMDD
    latest_file = None
    latest_date = None
    
    for log_file in log_files:
        # Extract the date from the filename using a regular expression
        match = re.search(r'(\d{4}_\d{2}_\d{2})', log_file)
        if match:
            file_date_str = match.group(0)
            file_date = datetime.strptime(file_date_str, date_format)
            
            # Compare to find the latest date
            if latest_date is None or file_date > latest_date:
                latest_date = file_date
                latest_file = log_file
        else:
            print("in find_latest_log_file function: No match found")
    
    return latest_file

# Load the CSV file into a DataFrame
csv_path = "/maps/epr26/tmf_pipe_out/project_memory_usage.csv"  # Replace with the path to your CSV file
if os.path.exists(csv_path):
    df = pd.read_csv(csv_path)
else:
    # Create an empty dataframe with expected columns if the CSV doesn't exist
    df = pd.DataFrame(columns=['proj_name', 'area_ha', 'full_acd', 'calculate_k_memory_GB', 'find_potential_matches_memory_GB', 
                               'build_m_table_memory_GB', 'find_pairs_memory_GB', 
                               'execution_time_seconds', 'execution_time_minutes'])

# Define the path to the directory containing the log files
project_directories = ["/maps/epr26/tmf_pipe_out_luc_t", "/maps/epr26/tmf_pipe_out_offset"]  # Replace with your directory name

new_project_paths = find_new_projects(project_directories)

# Iterate through new projects and extract information from their log files
for project_path in new_project_paths:
    # Retrieve directory and subfolder names
    project_path_split = os.path.split(project_path)
    project_directory = project_path_split[0]
    project_name = project_path_split[1]
    
    # Check if "additionality.csv" is present in the project folder
    additionality_file_path = os.path.join(project_path, "additionality.csv")
    if not os.path.exists(additionality_file_path):
        print(f"Skipping project '{project_name}' as 'additionality.csv' is not found.")
        memory_usage = {
            'calculate_k_memory_GB': np.nan,
            'find_potential_matches_memory_GB': np.nan,
            'build_m_table_memory_GB': np.nan,
            'find_pairs_memory_GB': np.nan
        }
        execution_time_seconds = np.nan
        execution_time_minutes = np.nan
    else:
        # Find the latest log file that matches the desired pattern
        log_file_path = find_latest_log_file(project_name, project_directory)
        
        if log_file_path and os.path.exists(log_file_path):
            # Parse the log file to retrieve memory usage and runtime
            memory_usage, execution_time_seconds = parse_log_file(log_file_path)
            
            # Calculate execution time in minutes
            execution_time_minutes = execution_time_seconds / 60 if execution_time_seconds else None
    
    # Retrieve project area value
    geojson_path = f"/maps/epr26/tmf-data/projects/{project_name}.geojson"
    with open(geojson_path) as f:
        gj = geojson.load(f)
        features = gj['features'][0]
        area_ha = area(features['geometry']) / 10000
    
    # Check availability of ACD for all LUC
    acd_path = os.path.join(project_path, "carbon-density.csv")
    if os.path.exists(acd_path):
        acd = pd.read_csv(acd_path)
        acd_class = {1, 2, 3, 4}
        full_acd = acd_class.issubset(acd.iloc[:, 0])
    else:
        full_acd = np.nan
    
    # Check if the data frame already contains a row for this project
    if project_name in df['proj_name'].values:
        # Fill in values to to the DataFrame if already present
        df.loc[df['proj_name'] == project_name, 'area_ha'] = area_ha
        df.loc[df['proj_name'] == project_name, 'full_acd'] = full_acd
        df.loc[df['proj_name'] == project_name, 'calculate_k_memory_GB'] = memory_usage['calculate_k_memory_GB']
        df.loc[df['proj_name'] == project_name, 'find_potential_matches_memory_GB'] = memory_usage['find_potential_matches_memory_GB']
        df.loc[df['proj_name'] == project_name, 'build_m_table_memory_GB'] = memory_usage['build_m_table_memory_GB']
        df.loc[df['proj_name'] == project_name, 'find_pairs_memory_GB'] = memory_usage['find_pairs_memory_GB']
        df.loc[df['proj_name'] == project_name, 'execution_time_seconds'] = execution_time_seconds
        df.loc[df['proj_name'] == project_name, 'execution_time_minutes'] = execution_time_minutes
    else:
        # Append new row to the DataFrame if not yet present
        new_row = pd.DataFrame([{
            'proj_name': project_name,
            'area_ha': area_ha,
            'full_acd': full_acd,
            'calculate_k_memory_GB': memory_usage['calculate_k_memory_GB'],
            'find_potential_matches_memory_GB': memory_usage['find_potential_matches_memory_GB'],
            'build_m_table_memory_GB': memory_usage['build_m_table_memory_GB'],
            'find_pairs_memory_GB': memory_usage['find_pairs_memory_GB'],
            'execution_time_seconds': execution_time_seconds,
            'execution_time_minutes': execution_time_minutes
        }])
        
        df = pd.concat([df, new_row], ignore_index = True)

df = df.sort_values('proj_name')

# Save the updated DataFrame back to the CSV file
df.to_csv(csv_path, index=False)