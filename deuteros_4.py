import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as mcolors

protein_length = 244
protein_name = '6s7e'
file_path = '250115_PCO4_Cluster_Table_FilteredPeptides.csv'

try:
    # Load the CSV file as a DataFrame
    df = pd.read_csv(file_path)
except FileNotFoundError:
    print(f"Error: File not found at {file_path}")
except pd.errors.EmptyDataError:
    print("Error: The CSV file is empty.")
except pd.errors.ParserError:
    print("Error: The CSV file could not be parsed.")

# Initialize an empty matrix to store the averaged values
unique_aa = range(1,protein_length)  # Adjust this based on the length of the protein
unique_exposure = sorted(df['Exposure'].unique())  # Ensure we have all exposure timepoints
# Initialize a dictionary to store the PyMOL scripts for each exposure timepoint
pymol_scripts = {}

# Find the global max absolute value of 'RU_diff_ave' across all exposure timepoints
min_ru = df['RU_diff_ave'].min()
max_ru = df['RU_diff_ave'].max()
max_abs = max(abs(min_ru), abs(max_ru))  # Symmetric max value

# Iterate over each exposure timepoint
for exposure in unique_exposure:
    # Create the PyMOL script for this exposure timepoint
    pymol_script = []

    # Get the average 'RU_diff_ave' values for each residue for the current exposure timepoint
    for aa_num in range(1, protein_length + 1):
        # Filter the dataframe for the current exposure and peptide containing the amino acid
        relevant_peptides = df[(df['Exposure'] == exposure)]
        relevant_peptides = relevant_peptides[(relevant_peptides['Start'] <= aa_num-1) & (relevant_peptides['End'] >= aa_num)] #removing the coverage for the first residue in the peptide
        
        # Filter out peptides that are significant (Protected or Deprotected)
        significant_peptides = relevant_peptides[relevant_peptides['Significance'].isin(['Protected', 'Deprotected'])]
        
        # If there are any significant peptides, average their 'RU_diff_ave' values
        if not significant_peptides.empty:
            avg_ru = significant_peptides['RU_diff_ave'].mean()
        else:
            avg_ru = 0  # Default value if no significant peptides exist
        # Normalize the average RU_diff_ave value symmetrically around 0
        normalized_avg_ru = avg_ru / max_abs if max_abs != 0 else 0        
        # Create a PyMOL command to set the B-factor for this residue
        pymol_script.append(f"alter resi {aa_num}, b = {normalized_avg_ru}")
    
        # Filter out peptides that are non significant
        non_significant_peptides = relevant_peptides[relevant_peptides['Significance'].isin(['Non Significant'])]

        # If there are any non_significant peptides, set their 'RU_diff_ave' value to 0
        if not non_significant_peptides.empty:
            avg_ru = 0  # Default value if no significant peptides exist
        # Normalize the average RU_diff_ave value symmetrically around 0
        normalized_avg_ru = avg_ru / max_abs if max_abs != 0 else 0
        # Create a PyMOL command to set the B-factor for this residue
        pymol_script.append(f"alter resi {aa_num}, b = {normalized_avg_ru}")
 
    # Combine the commands into a complete script for this exposure
    pymol_scripts[exposure] = "\n".join(pymol_script)

    # Optionally, write the PyMOL script to a file
    with open(f"colouring_script_exposure_{exposure}.pml", "w") as f:
        f.write("gray40, all\n")
        f.write(pymol_scripts[exposure])  # Write the B-factor adjustment commands
        f.write("\nrebuild\n")  # Rebuild to apply changes

    print(f"PyMOL script for exposure {exposure} written to colouring_script_exposure_{exposure}.pml")
