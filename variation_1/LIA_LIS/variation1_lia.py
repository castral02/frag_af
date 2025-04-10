import pandas as pd
import requests
from absl import app, flags, logging
import os
import re

# Define the output directory flag
flags.DEFINE_string('file', None, 'Output file name')
flags.DEFINE_string('excel', None, 'AlphaFold Metric Excel Sheet')
flags.DEFINE_string('uniprot', None, 'Uniprot ID of the tiled Protein ')
FLAGS = flags.FLAGS


def load_data(filepath):
    return pd.read_excel(filepath)

def extract_numeric(value):
    """
    Extract the numeric part from a given string.

    Args:
    - value (str): The string to extract the number from.

    Returns:
    - int: The extracted number or -1 if no number is found.
    """
    match = re.search(r"Fragment_(\d+)", value)
    if match:
        return int(match.group(1))  # Return as an integer
    return None  # Return None if no match is found

def filter_features(data, thresholds):
    for feature, threshold in thresholds.items():
        data = data[data[feature] >= threshold]
    return data

def get_protein_sequence(uniprot_id):
    """
    Fetch the protein sequence from UniProt for a given UniProt ID.

    Args:
    - uniprot_id (str): The UniProt ID of the protein.

    Returns:
    - str: The protein sequence.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        lines = response.text.split("\n")
        sequence = "".join(lines[1:])
        return sequence
    else:
        print(f"Failed to fetch protein sequence for UniProt ID {uniprot_id}. Status code: {response.status_code}")
        return None

def count_total_amino_acids(protein_sequence):
    """
    Count the total number of amino acids in the given protein sequence.

    Args:
    - protein_sequence (str): The protein sequence.

    Returns:
    - int: Total number of amino acids.
    """
    return len(protein_sequence)

def break_up_sequence(sequence_length, fragment_length, sliding_window):
    """
    Break up the protein sequence into fragments using a sliding window approach.

    Args:
    - sequence_length (int): Total length of the protein sequence.
    - fragment_length (int): Length of each fragment.
    - sliding_window (int): Step size of the sliding window.

    Returns:
    - list of tuples: Each tuple contains the start and end positions of a fragment.
    """
    fragments = []
    for start in range(1, sequence_length - fragment_length + 2, sliding_window):
        end = start + fragment_length - 1
        fragments.append((start, end))
    return fragments

def create_amino_acid_to_fragment_mapping(fragments):
    """
    Create a mapping of each amino acid to the fragments it belongs to.

    Args:
    - fragments (list of tuples): List of start and end positions of fragments.

    Returns:
    - dict: Mapping of amino acids to fragments.
    """
    amino_acid_to_fragment = {}
    for i, (start, end) in enumerate(fragments, start=1):
        for amino_acid in range(start, end + 1):
            if amino_acid not in amino_acid_to_fragment:
                amino_acid_to_fragment[amino_acid] = []
            amino_acid_to_fragment[amino_acid].append(i)
    return amino_acid_to_fragment

def replace_amino_acid_numbers_with_scores(amino_acid_to_fragment, amino_acid_scores):
    """
    Replace amino acid numbers with their corresponding scores from the fragments.

    Args:
    - amino_acid_to_fragment (dict): Mapping of amino acids to fragments.
    - amino_acid_scores (np.ndarray): Array of scores for each fragment.

    Returns:
    - list: List of scores for each amino acid.
    """
    replaced_array = []
    for amino_acid in sorted(amino_acid_to_fragment.keys()):
        fragment_numbers = amino_acid_to_fragment[amino_acid]
        scores = []
        for fragment_number in fragment_numbers:
            if fragment_number - 1 < len(amino_acid_scores):
                scores.append(amino_acid_scores[fragment_number - 1])
            else:
                print(f"Warning: Fragment number {fragment_number} is out of bounds for amino_acid_scores")
        replaced_array.append(scores)
    return replaced_array

def write_array_to_excel_with_adjusted_average(writer, array_data, sheet_name):
    """
    Write array data to an Excel sheet and adjust the average scores by subtracting the baseline.

    Args:
    - writer (pd.ExcelWriter): Excel writer object.
    - array_data (list): Array data to write.
    - sheet_name (str): Name of the Excel sheet.
    - baseline (float): Baseline value to subtract from the averages.
    """
    # Create a DataFrame from the array data
    df = pd.DataFrame(array_data)
    
    # Ensure all data is numeric, coerce non-numeric values to NaN
    df = df.apply(pd.to_numeric, errors='coerce')
    
    # Calculate row averages, skipping NaN values
    row_averages = df.mean(axis=1)
    
    # Add row averages as a new column
    df['Row Average'] = row_averages
    
    
    # Reorder columns to move the 'Row Average' and 'Data Minus Baseline' columns to the last positions
    columns = df.columns.tolist()
    columns = [col for col in columns if col not in ['Row Average']] + ['Row Average']
    df = df[columns]
    
    # Write the DataFrame to the specified sheet in the Excel file
    df.to_excel(writer, sheet_name=sheet_name, index=True)

def main (argv):
    excel_path = FLAGS.excel
    id = FLAGS.uniprot
    file = FLAGS.file

    data=load_data(excel_path)
    excel_df = load_data(excel_path)
    excel_df['numerical'] = excel_df['jobs'].apply(extract_numeric) # Extract the numeric part from the 'jobs' column and create a new 'numerical' column
    excel_df_sorted = excel_df.sort_values(by='numerical') # Sort the DataFrame by the 'numerical' column
    excel_df_sorted = excel_df_sorted.drop(columns=['numerical']) # Drop the 'numerical' column as it is no longer needed

    # Ensure the correct handling of headers and data rows
    excel_df_sorted.reset_index(drop=True, inplace=True)

    #filter data
    lia_threhsold = {
        "lis_score": 0.073,
        "average lia score": 1610
    }
    filter_df = filter_features(excel_df_sorted, lia_threhsold)

    #making a new column to grab the multiplier 
    excel_df_sorted['multiplier'] = excel_df_sorted['jobs'].apply(
    lambda x: -1 if x in filter_df['jobs'].values else 1)


    fragment_size = 60
    sliding_window = 10

    #grabbing protein
    protein_sequence = get_protein_sequence(id)
    if protein_sequence is None:
        return
    
    #multiplying columns by multiplier
    excel_df_sorted['iptm multiplied']=excel_df_sorted['iptm']*excel_df_sorted['multiplier']
    excel_df_sorted['mpdockq multiplied']=excel_df_sorted['mpDockQ/pDockQ']*excel_df_sorted['multiplier']

    #grabbing specific columns
    iptm_m = excel_df_sorted['iptm multiplied']
    mpdock_m =excel_df_sorted['mpdockq multiplied']

    total_amino_acids = count_total_amino_acids(protein_sequence)
    
    fragments = break_up_sequence(total_amino_acids, fragment_size, sliding_window)
    
    amino_acid_to_fragment = create_amino_acid_to_fragment_mapping(fragments)

    replaced_array_iptm = replace_amino_acid_numbers_with_scores(amino_acid_to_fragment, iptm_m)
    replaced_array_mpdock = replace_amino_acid_numbers_with_scores(amino_acid_to_fragment, mpdock_m)
    
    path = '/Users/castroverdeac/Desktop/codes_for_AF/actual_frag_af/data/variation1_prediction_lia/'
    output_file=os.path.join(path, f"{file}.xlsx")
    with pd.ExcelWriter(output_file) as writer:
        write_array_to_excel_with_adjusted_average(writer, replaced_array_iptm, 'iPTM Scores')
        write_array_to_excel_with_adjusted_average(writer, replaced_array_mpdock, 'mpDockQ Scores')
    
    print(f"Array data has been written to {output_file}")
   
if __name__ == '__main__':
    app.run(main)
