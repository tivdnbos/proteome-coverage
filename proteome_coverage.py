import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Define the function to calculate proteome coverage
def calculate_coverage(peptide_df, fasta_file_path):
    # Load the FASTA file
    with open(fasta_file_path, 'r') as f:
        fasta_lines = f.readlines()

    # Parse the FASTA file to obtain the protein sequences
    protein_sequences = {}
    current_protein = ''
    for line in fasta_lines:
        if line.startswith('>'):
            current_protein = line[1:].strip()
            protein_sequences[current_protein] = ''
        else:
            protein_sequences[current_protein] += line.strip()

    # Calculate the length of each protein
    protein_lengths = {}
    for protein in protein_sequences:
        protein_lengths[protein] = len(protein_sequences[protein])

    # Calculate the number of peptides that cover each protein
    protein_coverage = {}
    for protein in protein_sequences:
        protein_coverage[protein] = []
    for index, row in peptide_df.iterrows():
        sequence = row['peptide']
        for protein in protein_sequences:
            if sequence in protein_sequences[protein]:
                start = protein_sequences[protein].find(sequence)
                covered_positions = np.arange(start, start + len(sequence), 1).tolist()
                protein_coverage[protein].extend(covered_positions)
    print(protein_coverage)

    # Calculate the percent coverage of each protein
    protein_percent_coverage = {}
    for protein in protein_sequences:
        protein_percent_coverage[protein] = (len(set(protein_coverage[protein])) / protein_lengths[protein]) * 100

    # Convert the dictionary of percent coverage values to a DataFrame
    coverage_df = pd.DataFrame.from_dict(protein_percent_coverage, orient='index', columns=['coverage'])
    print(coverage_df)
    return coverage_df

# Define the function to plot the proteome coverage distribution
def plot_coverage(coverage_df, plot_title, plot_filename):
    # Create the plot
    fig, ax = plt.subplots()
    sns.histplot(data=coverage_df, x='coverage', binwidth=5, binrange=(0, 100), kde=True, ax=ax)

    # Set the plot title and labels
    ax.set_title(plot_title)
    ax.set_xlabel('Percent coverage')
    ax.set_ylabel('Number of proteins')

    # Save the plot to a file
    fig.savefig(plot_filename)

# Load the peptide data
peptide_df = pd.read_csv('/public/compomics2/Caro_BactPTMs/PRIDE_MS2Rescore_results/K12_ref/PXD000498/A14-07066_searchengine_ms2pip_rt_features_protformat.pout', sep = "\t")
print(peptide_df)
# Define the path to the FASTA file
fasta_file_path = '/home/compomics/Documents/ionbot_workdir/ms2rescore_reprocessing/reference_dbs/uniprot_UP000000625_17-11-22_K12_with_crap.fasta'

# Calculate the proteome coverage
coverage_df = calculate_coverage(peptide_df, fasta_file_path)

# Plot the proteome coverage distribution
plot_title = 'Proteome coverage distribution'
plot_filename = 'proteome_coverage_distribution.png'
plot_coverage(coverage_df, plot_title, plot_filename)