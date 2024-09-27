import os
import re
import pandas as pd
from glob import glob
from datetime import datetime

# 1. Load the CSV file into a DataFrame
csv_path = "/maps/epr26/tmf_pipe_out/project_memory_usage.csv"  # Replace with the path to your CSV file
if os.path.exists(csv_path):
    df = pd.read_csv(csv_path)
else:
    # Create an empty dataframe with expected columns if the CSV doesn't exist
    df = pd.DataFrame(columns=['proj_name', 'calculate_k_memory_GB', 'find_potential_matches_memory_GB', 
                               'build_m_table_memory_GB', 'find_pairs_memory_GB', 
                               'execution_time_seconds', 'execution_time_minutes'])

# Define the path to the directory containing the log files
log_directory = "/maps/epr26/tmf_pipe_out_luc_t"  # Replace with your directory name

# 1. Get list of existing project names in the DataFrame
existing_projects = set(df['proj_name'].tolist())

# 1. Get list of project names based on subfolders in log_directory
project_subfolders = set(os.listdir(log_directory))
project_subfolders = {entry.name for entry in os.scandir(log_directory) if entry.is_dir()}

# 1.1 Filter out subfolders that are of format "as[0-9][0-9]" or start with "rescaled"
filtered_subfolders = set()
for folder in project_subfolders:
    if re.match(r"as\d{2}$", folder):  # Exclude "as[0-9][0-9]"
        continue
    if folder.startswith(("slopes", "rescaled", "srtm")):  # Exclude "slopes" and "rescaled"
        continue
    filtered_subfolders.add(folder)

# Find new projects to process (those not already in the CSV)
new_projects = filtered_subfolders - existing_projects

# 2. Function to parse memory usage and runtime from log files
def parse_log_file(log_file_path):
    memory_usage = {
        'calculate_k_memory_GB': None,
        'find_potential_matches_memory_GB': None,
        'build_m_table_memory_GB': None,
        'find_pairs_memory_GB': None
    }
    execution_time_seconds = None
    
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
    print(log_files)
    
    # Check if any files match the pattern
    if not log_files:
        return None

    # Extract date from each file name and find the latest one
    date_format = "%Y_%m_%d"  # Assuming date in filename is in format YYYYMMDD
    latest_file = None
    latest_date = None
    
    for log_file in log_files:
        # Extract the date from the filename using a regular expression
        match = re.search(r'(\d{4}_\d{2}_\d{2})', log_file)
        if match:
            print(f"Match found: {match.group(0)}")
            #@@@it still doesn't find a match
            file_date_str = match.group(1)
            file_date = datetime.strptime(file_date_str, date_format)
            
            # Compare to find the latest date
            if latest_date is None or file_date > latest_date:
                latest_date = file_date
                latest_file = log_file
    
    return latest_file


# 5. Iterate through new projects and extract information from their log files
for project_name in new_projects:
    # Construct the path to the project folder
    project_folder_path = os.path.join(log_directory, project_name)
    
    # Check if "additionality.csv" is present in the project folder
    additionality_file_path = os.path.join(project_folder_path, "additionality.csv")
    if not os.path.exists(additionality_file_path):
        print(f"Skipping project '{project_name}' as 'additionality.csv' is not found.")
        continue  # Skip this project as it may still be running
    
    # Find the latest log file that matches the desired pattern
    log_file_path = find_latest_log_file(project_name, log_directory)
    print("log_file_path:"+str(log_file_path))
    
    if log_file_path and os.path.exists(log_file_path):
        # 2. Parse the log file to retrieve memory usage and runtime
        memory_usage, execution_time_seconds = parse_log_file(log_file_path)
        print(f"execution_time_seconds: {execution_time_seconds}")
        
        # 4. Calculate execution time in minutes
        execution_time_minutes = execution_time_seconds / 60 if execution_time_seconds else None
        
        # 5. Append new information to the DataFrame
        df = df.append({
            'proj_name': project_name,
            'calculate_k_memory_GB': memory_usage['calculate_k_memory_GB'],
            'find_potential_matches_memory_GB': memory_usage['find_potential_matches_memory_GB'],
            'build_m_table_memory_GB': memory_usage['build_m_table_memory_GB'],
            'find_pairs_memory_GB': memory_usage['find_pairs_memory_GB'],
            'execution_time_seconds': execution_time_seconds,
            'execution_time_minutes': execution_time_minutes
        }, ignore_index=True)

# Save the updated DataFrame back to the CSV file
df.to_csv(csv_path, index=False)