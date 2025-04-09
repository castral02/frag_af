#in this variation, I am going to have all three in this code. But also, the we would find the average mpDockQ score and multiple by mulitplier that was found
import pandas as pd
import requests
from absl import app, flags, logging
import os
import re
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import Birch
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Define the output directory flag
flags.DEFINE_string('file', None, 'Output file name')
flags.DEFINE_string('excel', None, 'AlphaFold Metric Excel Sheet')
flags.DEFINE_string('uniprot', None, 'Uniprot ID of the tiled Protein ')
FLAGS = flags.FLAGS


def load_data(filepath):
    return pd.read_excel(filepath)

def filter_features(data, thresholds):
    for feature, threshold in thresholds.items():
        data = data[data[feature] >= threshold]
    return data
def extract_numeric(job_name):
    match = re.search(r"Fragment_(\d+)", job_name)
    if match:
        return int(match.group(1))  # Return as an integer
    return None  # Return None if no match is found
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


def write_array_to_excel_with_adjusted_average(array_data, baseline_df):
    """
    Write array data to an Excel sheet and adjust the average scores by multiplying the baseline DataFrame.

    Args:
    - array_data (list): Array data to write.
    - baseline_df (pd.DataFrame): Baseline DataFrame to multiply from the data.
    """
    # Create a DataFrame from the array data
    df = pd.DataFrame(array_data)
    
    # Ensure all data is numeric, coerce non-numeric values to NaN
    df = df.apply(pd.to_numeric, errors='coerce')
    
    # Align baseline DataFrame with the main DataFrame by index and columns
    baseline_df = baseline_df.reindex_like(df)
    
    # Calculate row averages, skipping NaN values
    row_averages = df.mean(axis=1)
    
    # Add row averages as a new column
    df['Row Average'] = row_averages
    df['Multipler'] = baseline_df [0]
    
    # Subtract the baseline from the data and calculate the new column
    df['Times Mulitplier'] = row_averages*baseline_df[0]
    
    # Reorder columns to move the 'Row Average' and 'Data Minus Baseline' columns to the last positions
    columns = df.columns.tolist()
    df = df[columns]

    return df 

def main (argv):
    excel_path = FLAGS.excel
    id = FLAGS.uniprot
    file = FLAGS.file

    #step 1: load data
    try:
        data_df = load_data(excel_path)
        logging.info(f"Successfully loaded Excel file from {excel_path}.")
    except FileNotFoundError:
        logging.error(f"Excel file not found at {excel_path}.")
        return
    except Exception as e:
        logging.error(f"Error loading Excel file: {e}")
        return
    
    #step 2: cleaning dataframe
    data_df['numerical'] = data_df['jobs'].apply(extract_numeric)     # Extract the numeric part from the 'jobs' column and create a new 'numerical' column
    data_df_sorted = data_df.sort_values(by='numerical')    # Sort the DataFrame by the 'numerical' column
    data_df_sorted = data_df_sorted.drop(columns=['numerical'])    # Drop the 'numerical' column as it is no longer needed
    data_df_sorted.reset_index(drop=True, inplace=True)    # Ensure the correct handling of headers and data rows

    sorted_mpdock = data_df_sorted.copy(deep=True) #this is what you are going to put in the excel
    sorted_lia =data_df_sorted.copy(deep=True) 
    #step 3: filter
    lia_threhsold = {
        "mpDockQ/pDockQ": 0.175,
    }
    mpdock_threhshold = {
       "iptm": 0.5
    }
    lia_df = filter_features(sorted_lia, lia_threhsold)
    mpdock_df = filter_features(sorted_mpdock, mpdock_threhshold)
    if lia_df.empty:
        logging.warning("No data passed the filtering LIA/LIS threshold.")
    else:
        logging.info(f"Filtered data contains {len(lia_df)} rows.")
    
    if mpdock_df.empty:
        logging.warning("No data passed the filtering mpDockQ/contact threshold.")
    else:
        logging.info(f"Filtered data contains {len(mpdock_df)} rows.")

    #step 4: normalizing data
    columns_mpdock = ['average pae score', 'Hydrophobhic', 'int_solv_en', 'Num_intf_residues', 'int_area']
    columns_lis = ['pi_score', 'mpDockQ/pDockQ', 'Polar', 'Num_intf_residues', 'contact_pairs']

    missing_cols = [col for col in columns_mpdock if col not in mpdock_df.columns]
    if missing_cols:
        logging.error(f"Missing columns for normalization: {missing_cols}")
        return
    else:
        logging.info("All required columns for normalization are present.")
    
    missing_cols = [col for col in columns_lis if col not in lia_df.columns]
    if missing_cols:
        logging.error(f"Missing columns for normalization: {missing_cols}")
        return
    else:
        logging.info("All required columns for normalization are present.")

    "handling empty filtered data"
    if mpdock_df.empty:
        normdock_df = pd.DataFrame(columns=mpdock_df.columns)
        logging.info("Created an empty DataFrame for normalization since no data passed filtering.")
    else:
        normdock_df = mpdock_df.copy()
        scaler = StandardScaler()
        
        # Normalize specified columns
        normdock_df[columns_mpdock] = scaler.fit_transform(mpdock_df[columns_mpdock])
        logging.info("Data normalization completed successfully.")
        
        # Rename columns
        rename_dict = {
            'average pae score': 'norm_pae',
            'Hydrophobhic':'norm_hydro',
            'int_solv_en': 'norm_int',
            'Num_intf_residues': 'norm_num',
            'int_area':'norm_int_a'
        }
        normdock_df.rename(columns=rename_dict, inplace=True)

    "handling empty filtered data"
    if lia_df.empty:
        normlia_df = pd.DataFrame(columns=lia_df.columns)
        logging.info("Created an empty DataFrame for normalization since no data passed filtering.")
    else:
        normlia_df = lia_df.copy()
        scaler = StandardScaler()
        
        # Normalize specified columns
        normlia_df[columns_lis] = scaler.fit_transform(lia_df[columns_lis])
        logging.info("Data normalization completed successfully.")
        
        # Rename columns
        rename_dict = {
            'pi_score':'norm_pi',
            'mpDockQ/pDockQ':'norm_mpdock',
            'Polar': 'norm_polar',
            'Num_intf_residues': 'norm_num',
            'contact_pairs':'norm_contact'
        }
        normlia_df.rename(columns=rename_dict, inplace=True)

    #Step 4: Clustering
    if not normdock_df.empty:
       birch = Birch(branching_factor=50, threshold=0.3, n_clusters=4)
       birch.fit(normdock_df[['norm_pae', 'norm_num', 'norm_int', 'norm_int_a', 'norm_hydro']])
       labels = birch.labels_
       normdock_df['Cluster'] = labels
       logging.info(f"Clustering completed with {len(set(labels))} clusters identified.")
    else:
       logging.warning("Skipping clustering due to empty normalized data.")

    if not normlia_df.empty:
       birch = Birch(branching_factor=50, threshold=0.3, n_clusters=4)
       birch.fit(normlia_df[['norm_mpdock', 'norm_pi', 'norm_num', 'norm_contact', 'norm_polar']])       
       labels = birch.labels_
       normlia_df['Cluster'] = labels
       logging.info(f"Clustering completed with {len(set(labels))} clusters identified.")
    else:
       logging.warning("Skipping clustering due to empty normalized data.")

    #Step 5: grabbing protein sequence and finding average mpDockQ
    sequence = get_protein_sequence(id)
    if sequence is None:
       logging.error(f"Failed to fetch protein sequence for UniProt ID {id}.")
       return
    else:
       logging.info(f"Protein sequence fetched successfully with length {len(sequence)}.")
    
    fragment_size = 60
    sliding_window = 10
    mpdock_column=data_df_sorted['mpDockQ/pDockQ']
    total_amino_acids = count_total_amino_acids(sequence)
    fragments = break_up_sequence(total_amino_acids, fragment_size, sliding_window)
    amino_acid_to_fragment = create_amino_acid_to_fragment_mapping(fragments)
    average_mpDockQ = replace_amino_acid_numbers_with_scores(amino_acid_to_fragment, mpdock_column)
    average_mpDockQ =pd.DataFrame(average_mpDockQ)
    lis_copy = average_mpDockQ.copy(deep=True) #this is where you append the new materials
    mpdock_copy = average_mpDockQ.copy(deep=True)

    #Step 6: Mapping to the original DataFrame
    if 'Cluster' not in normdock_df.columns: 
        normdock_df['Predicted_Priority']='not_pass'
        job_pred_priority_mpdock = dict(zip(normdock_df['jobs'], normdock_df['Predicted_Priority']))
        sorted_mpdock['Predicted_Priority'] = sorted_mpdock['jobs'].map(job_pred_priority_mpdock).fillna('not_pass')
        print("The 'cluster' column is missing.")
    else:
        try:
            label_mapping = {0: 'not_pass', 1: 'Low', 2: 'Medium', 3: 'High'}
            labels_mpdock = normdock_df['Cluster']  # or any other way to define labels
            if 'jobs' not in normdock_df.columns or 'jobs' not in sorted_mpdock.columns:
                logging.error("'jobs' column missing in normalized_df or data_df_sorted.")
                raise KeyError("'jobs' column is required.")
            try:
                normdock_df['Predicted_Priority'] = [label_mapping[pred] if pred in label_mapping else 'Unknown' for pred in labels_mpdock]
            except KeyError as e:
                logging.error(f"Mapping error: {e}")
                raise
            job_pred_priority_mpdock = dict(zip(normdock_df['jobs'], normdock_df['Predicted_Priority']))
            sorted_mpdock['Predicted_Priority'] = sorted_mpdock['jobs'].map(job_pred_priority_mpdock).fillna('not_pass')
            logging.info("Successfully mapped predicted priorities to the original dataframe.")
        except KeyError as e:
            logging.error(f"Mapping error: {e}")
            raise
        except Exception as e:
            logging.error(f"Unexpected error during mapping: {e}")
            raise
    
    # Visualization
    priority_map = {
        'not_pass': 1,
        'low': -1,
        'medium': -2,
        'high': -3
    }

    # Match fragments to tiled sequences
    job_titles_mpdock = sorted_mpdock['jobs']  # Extract jobs column
    fragment_numbers_mpdock = job_titles_mpdock.str.extract(r'Fragment_(\d+)', expand=False)  # Extract fragment number
    fragment_numbers_mpdock = fragment_numbers_mpdock.astype(int) - 1  # Convert to zero-based index
    
    # Add corresponding tiled sequences to the processed_data DataFrame
    sorted_mpdock['tiled_sequence'] = fragment_numbers_mpdock.apply(
        lambda idx: sequence[idx*10:idx*10+60] if idx*10 + 60 <= len(sequence) else None
    )
    sequence_length = len(sequence)
    amino_acid_scores_mpdock = np.zeros(sequence_length)

    for _, fragment in sorted_mpdock.iterrows():
        fragment_sequence = fragment['tiled_sequence']
        if fragment_sequence:
            start = sequence.find(fragment_sequence)  # Find the start position of the fragment in the sequence
            end = start + len(fragment_sequence)

            # Map the priority to the sequence
            score_mpdock = priority_map.get(fragment['Predicted_Priority'].lower(), 0)  # Default to 0 if no match
            amino_acid_scores_mpdock[start:end] += score_mpdock / (end - start)  # Distribute the score evenly
    "lis/lia"
    if 'Cluster' not in normlia_df.columns: 
        normlia_df['Predicted_Priority']='not_pass'
        job_pred_priority_lia = dict(zip(normlia_df['jobs'], normlia_df['Predicted_Priority']))
        sorted_lia['Predicted_Priority'] = sorted_lia['jobs'].map(job_pred_priority_lia).fillna('not_pass')
        print("The 'cluster' column is missing.")
    else:
        try:
            label_mapping = {0: 'not_pass', 1: 'Low', 2: 'Medium', 3: 'High'}
            labels_lia = normlia_df['Cluster']  # or any other way to define labels
            if 'jobs' not in normlia_df.columns or 'jobs' not in sorted_lia.columns:
                logging.error("'jobs' column missing in normalized_df or data_df_sorted.")
                raise KeyError("'jobs' column is required.")
            try:
                normlia_df['Predicted_Priority'] = [label_mapping[pred] if pred in label_mapping else 'Unknown' for pred in labels_lia]
            except KeyError as e:
                logging.error(f"Mapping error: {e}")
                raise
            job_pred_priority_lia = dict(zip(normlia_df['jobs'], normlia_df['Predicted_Priority']))
            sorted_lia['Predicted_Priority'] = sorted_lia['jobs'].map(job_pred_priority_lia).fillna('not_pass')
            logging.info("Successfully mapped predicted priorities to the original dataframe.")
        except KeyError as e:
            logging.error(f"Mapping error: {e}")
            raise
        except Exception as e:
            logging.error(f"Unexpected error during mapping: {e}")
            raise
    
    # Visualization
    priority_map = {
        'not_pass': 1,
        'low': -1,
        'medium': -2,
        'high': -3
    }

    # Match fragments to tiled sequences
    job_titles_lia = sorted_lia['jobs']  # Extract jobs column
    fragment_numbers_lia = job_titles_lia.str.extract(r'Fragment_(\d+)', expand=False)  # Extract fragment number
    fragment_numbers_lia = fragment_numbers_lia.astype(int) - 1  # Convert to zero-based index
    
    # Add corresponding tiled sequences to the processed_data DataFrame
    sorted_lia['tiled_sequence'] = fragment_numbers_lia.apply(
        lambda idx: sequence[idx*10:idx*10+60] if idx*10 + 60 <= len(sequence) else None
    )
    sequence_length = len(sequence)
    amino_acid_scores_lia = np.zeros(sequence_length)

    for _, fragment in sorted_lia.iterrows():
        fragment_sequence = fragment['tiled_sequence']
        if fragment_sequence:
            start = sequence.find(fragment_sequence)  # Find the start position of the fragment in the sequence
            end = start + len(fragment_sequence)

            # Map the priority to the sequence
            score_lia = priority_map.get(fragment['Predicted_Priority'].lower(), 0)  # Default to 0 if no match
            amino_acid_scores_lia[start:end] += score_lia / (end - start)  # Distribute the score evenly
    amino_acid_scores_mpdock = pd.DataFrame(amino_acid_scores_mpdock)
    amino_acid_scores_lia = pd.DataFrame(amino_acid_scores_lia)

    #Step 7: multiplying number to mpDockQ 
    lis_excel = write_array_to_excel_with_adjusted_average(lis_copy, amino_acid_scores_lia)
    mpdock_excel = write_array_to_excel_with_adjusted_average(mpdock_copy, amino_acid_scores_mpdock)

      # Step 11: Excel output
    try:
        path_output = '/Users/castroverdeac/Desktop/codes_for_AF/actual_frag_af/data/variation2_predictions_individual'
        output_file = os.path.join(path_output, f"{file}.xlsx")
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            data_df_sorted.to_excel(writer, sheet_name='predictions', index=False)
            normdock_df.to_excel(writer, sheet_name='iptm', index=False)
            normlia_df.to_excel(writer, sheet_name='mpdock', index=False)
            lis_excel.to_excel(writer, sheet_name='mpdock scores', index=True)
            mpdock_excel.to_excel(writer, sheet_name='iptm scores', index=True)
        logging.info(f"Excel file successfully saved at {output_file}.")
    except Exception as e:
        logging.error(f"Failed to save Excel file: {e}")
    
if __name__ == '__main__':
    app.run(main)