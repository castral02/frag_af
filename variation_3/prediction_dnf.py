import joblib
import pandas as pd
import requests
from absl import flags, app, logging
import numpy as np
from matplotlib.colors import Normalize
import re
import os
from sklearn.preprocessing import MinMaxScaler


#-------------- Defining Flags --------------#
flags.DEFINE_string('file', None, 'Output file name')
flags.DEFINE_string('excel', None, 'AlphaFold Metric Excel Sheet')
flags.DEFINE_string('uniprot', None, 'Uniprot ID of the tiled Protein ')
FLAGS = flags.FLAGS

#-------------- Definitions --------------#
def load_data(filepath):
    return pd.read_excel(filepath, sheet_name='brafT_brafF')

def filter_features(data, thresholds):
    for feature, threshold in thresholds.items():
        data = data[data[feature] >= threshold]
    numeric_cols = ['Polar',
                         'contact_pairs','sc',
                         'Hydrophobhic',
                         'average pae score', 'mpDockQ/pDockQ']
    data = data.dropna(subset=numeric_cols)
    return data

def normalize_drop(data):
    data.columns = data.columns.str.strip()
    numeric_cols = ['Polar',
                         'contact_pairs','sc',
                         'Hydrophobhic',
                         'average pae score', 'mpDockQ/pDockQ']
    data[numeric_cols] = data.groupby('tiled protein')[numeric_cols].transform(
        lambda x: MinMaxScaler().fit_transform(x.values.reshape(-1, 1)).flatten()
    )
    return data

def extract_numeric(job_name):
    match = re.search(r"Fragment_(\d+)", job_name)
    if match:
        return int(match.group(1))  # Return as an integer
    return None  # Return None if no match is found

def make_predictions(model, X_new):
    """Use the loaded model to make predictions."""
    predictions = model.predict(X_new)
    return predictions

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
    df = pd.DataFrame(array_data)
    df = df.apply(pd.to_numeric, errors='coerce')
    baseline_df = pd.DataFrame(baseline_df)
    baseline_df = baseline_df.apply(pd.to_numeric, errors='coerce')
    baseline_df = baseline_df.reindex_like(df)
    row_averages = df.mean(axis=1)
    df['Row Average'] = row_averages
    df['Multipler'] = baseline_df [0]
    df['Times Mulitplier'] = row_averages*baseline_df[0]
    columns = df.columns.tolist()
    df = df[columns]
    return df 

def main(argv):
    excel_path = FLAGS.excel
    id = FLAGS.uniprot
    file = FLAGS.file

#-------------- S1: Load Data --------------#
    try:
        data_df = load_data(excel_path)
        logging.info(f"Successfully loaded Excel file from {excel_path}.")
    except FileNotFoundError:
        logging.error(f"Excel file not found at {excel_path}.")
        return
    except Exception as e:
        logging.error(f"Error loading Excel file: {e}")
        return
    
#-------------- S2: Load Model --------------#
    try:
        mpdock_path = '/Users/castroverdeac/Desktop/codes_for_AF/actual_frag_af/model_variations/variation_3/XGBoost/model/model_xgboost_model.pkl'
        mpdock_model= joblib.load(mpdock_path)
        print(f"Loaded model from {mpdock_path}")
    except FileNotFoundError:
        logging.error(f"Model not found at {mpdock_path}.")
        return
    except Exception as e:
        logging.error(f"Error loading model: {e}")
        return

#-------------- S3: Sorting DataFrame --------------#
    data_df['numerical'] = data_df['jobs'].apply(extract_numeric)     # Extract the numeric part from the 'jobs' column and create a new 'numerical' column    
    sorted_mpdock = data_df.copy(deep=True) #this is what you are going to put in the excel

#-------------- S4: Data Proceessing --------------#
    mpdock_threhshold = {'iptm': 0.5}
    mpdock_df = filter_features(sorted_mpdock, mpdock_threhshold)
    normalize_features = normalize_drop(mpdock_df)
    cluster_features_mpdockq = ['Polar',
                         'contact_pairs','sc',
                         'Hydrophobhic',
                         'average pae score', 'mpDockQ/pDockQ']
    processed_mpdock= normalize_features[cluster_features_mpdockq]

#-------------- S5: Making Prediction --------------#
    if mpdock_df.empty:
        logging.warning("The 'mpdock_df' DataFrame is empty. No predictions will be made for 'mpdock_model'.")
        mpdock_predictions=[]
    else:
        mpdock_predictions = make_predictions(mpdock_model, processed_mpdock)
        logging.info(f"Predictions made for 'mpdock_model': {len(mpdock_predictions)} entries.")
    
    label_mapping = {0: 'not_pass', 1: 'Low', 2: 'Medium', 3: 'High'}

#-------------- S6: Mapping Predicion to DataFrame  --------------#
    if len(mpdock_predictions) == 0:
        mpdock_df['Predicted_Priority']='not_pass'
        job_pred_priority_mpdock = dict(zip(mpdock_df['jobs'], mpdock_df['Predicted_Priority']))
        sorted_mpdock['Predicted_Priority'] = sorted_mpdock['jobs'].map(job_pred_priority_mpdock).fillna('not_pass')
    else:
        mpdock_df['Predicted_Priority'] = [label_mapping[pred] for pred in mpdock_predictions]

#-------------- S7: Grabbing Protein sequence and Finding Average mpDockQ  --------------#
    sequence = get_protein_sequence(id)
    if sequence is None:
       logging.error(f"Failed to fetch protein sequence for UniProt ID {id}.")
       return
    else:
       logging.info(f"Protein sequence fetched successfully with length {len(sequence)}.")
    
    fragment_size = 60
    sliding_window = 10
    mpdock_column=data_df['mpDockQ/pDockQ']
    total_amino_acids = count_total_amino_acids(sequence)
    fragments = break_up_sequence(total_amino_acids, fragment_size, sliding_window)
    amino_acid_to_fragment = create_amino_acid_to_fragment_mapping(fragments)
    average_mpDockQ = replace_amino_acid_numbers_with_scores(amino_acid_to_fragment, mpdock_column)
    average_mpDockQ =pd.DataFrame(average_mpDockQ)

    mpdock_copy = average_mpDockQ.copy(deep=True)

#-------------- S8:Mapping to Original DataFrame --------------#
    try:
        label_mapping = {0: 'not_pass', 1: 'Low', 2: 'Medium', 3: 'High'}
        if 'jobs' not in mpdock_df.columns or 'jobs' not in sorted_mpdock.columns:
            logging.error("'jobs' column missing in normalized_df or data_df_sorted.")
            raise KeyError("'jobs' column is required.")
        job_pred_priority_mpdock = dict(zip(mpdock_df['jobs'], mpdock_df['Predicted_Priority']))
        sorted_mpdock['Predicted_Priority'] = sorted_mpdock['jobs'].map(job_pred_priority_mpdock).fillna('not_pass')
        logging.info("Successfully mapped predicted priorities to the original dataframe.")
    except KeyError as e:
        logging.error(f"Mapping error: {e}")
        raise
    except Exception as e:
        logging.error(f"Unexpected error during mapping: {e}")
        raise

#-------------- S9: Visualization  --------------#
    priority_map = {
        'not_pass': 1,
        'low': -1,
        'medium': -2,
        'high': -3
    }

    job_titles_mpdock = sorted_mpdock['jobs']  # Extract jobs column
    fragment_numbers_mpdock = job_titles_mpdock.str.extract(r'Fragment_(\d+)', expand=False)  # Extract fragment number
    fragment_numbers_mpdock = fragment_numbers_mpdock.astype(int) - 1  # Convert to zero-based index
    
    sorted_mpdock['tiled_sequence'] = fragment_numbers_mpdock.apply(
        lambda idx: sequence[idx*10:idx*10+60] if idx*10 + 60 <= len(sequence) else None
    )
    sequence_length = len(sequence)
    amino_acid_scores_mpdock = np.zeros(sequence_length)

    for _, fragment in sorted_mpdock.iterrows():
        fragment_sequence = fragment['tiled_sequence']
        if fragment_sequence:
            start = sequence.find(str(fragment_sequence))  # Ensure fragment_sequence is treated as a string
            end = start + len(fragment_sequence)

            # Map the priority to the sequence
            score_mpdock = priority_map.get(fragment['Predicted_Priority'].lower(), 0)  # Default to 0 if no match
            amino_acid_scores_mpdock[start:end] += score_mpdock / (end - start)  # Distribute the score evenly

#-------------- S10: multiplying number to mpDockQ  --------------#
    mpdock_excel = write_array_to_excel_with_adjusted_average(mpdock_copy, amino_acid_scores_mpdock)

#-------------- S11:Excel Output  --------------#
    try:
        path_output = '/Users/castroverdeac/Desktop/codes_for_AF/actual_frag_af/model_variations/variation_3/XGBoost/log'
        output_file = os.path.join(path_output, f"{file}.xlsx")
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            data_df.to_excel(writer, sheet_name='predictions', index=False)
            mpdock_df.to_excel(writer, sheet_name='normalized mpDockQ Contact data', index=False)
            mpdock_excel.to_excel(writer, sheet_name='mpDockQ scores', index=True)
        logging.info(f"Excel file successfully saved at {output_file}.")
    except Exception as e:
        logging.error(f"Failed to save Excel file: {e}")

if __name__ == '__main__':
    app.run(main)
