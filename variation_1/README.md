# Variation 1: Correlating AlphaFold Metrics with DNFs

As initial approaches, we correlated commonly used AlphaFold Metrics to dominant negative fragments. In this repository, we explored 3 metrics: Local Interaction Area and Local Interaction Score, mpDockQ, and ipTM.

## mpDockQ & ipTM
We explored two popular metrics: interface predicted templating modelling ([ipTM](https://academic.oup.com/bioinformatics/article/26/7/889/213219?login=true)) and [mpDockQ](https://www.nature.com/articles/s41467-022-33729-4).

Thresholds were determined from previous research:

- mpDockQ: ≥ 0.173 (1, 3)
- ipTM: ≥ 0.5 (2, 5)

**Workflow Image**

![Workflow](../images/workflow/frag_af_variation_1.heic)

1. AlphaPulldown Workflow
2. Extracting AlphaFold Metrics

mpDockQ/ipTM
```python
iptm_column = excel_df_sorted['iptm'].astype(float).values  # Assuming 'iptm' is the column name
mpdock_column = excel_df_sorted['mpDockQ/pDockQ'].astype(float).values  # Assuming 'mpDockQ/pDockQ' is the column name
```

4. Mapping Metrics onto individual amino acids

mpDockQ/ipTM
```python
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
```

5. Averaging and Subtracting Baseline

mpDockQ/ipTM
```python
def write_array_to_excel_with_adjusted_average(writer, array_data, sheet_name, baseline):
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
    
    # Add a new column with the average minus the baseline
    df['Data Minus Baseline'] =  baseline - df['Row Average']
    
    # Reorder columns to move the 'Row Average' and 'Data Minus Baseline' columns to the last positions
    columns = df.columns.tolist()
    columns = [col for col in columns if col not in ['Row Average', 'Data Minus Baseline']] + ['Row Average', 'Data Minus Baseline']
    df = df[columns]
    
    # Write the DataFrame to the specified sheet in the Excel file
    df.to_excel(writer, sheet_name=sheet_name, index=True)
```

8. Inferring Dominant Negative Fragments (DNFs)

### Dependencies to download
```bash
pip install pandas numpy absl-py biopython openpyxl
```
You can also [click here](mpdockq_iptm/variation_1_mpdockq_iptm.yml) to download the conda environment.

### How to run: 
This is the [code](mpdockq_iptm/variation_1_iptm_mpdock.py) to run.
```bash
python3 variation_1.py -uniprot=uniprot_id -file=excel_output_name -excel=/path/to/AlphaPulldown/outputs
```

Here is an example of the [output](../pipeline/example/flt3_iptm_mpdockq_v1.xlsx)

## LIA & LIS

Local Interaction Area (LIA) and Local Interaction Score (LIS) is a newly developed metric from the Perrimon Lab at Harvard. The scoring system is derived off of the AlphaFold Metric Predicted Aligned Error (PAE) to discover highly interactive areas when the strucutre is small and felxible (4). 

We used the average LIA/LIS scores of each fragment and the threshold's from the original paper-- *Enhanced Protein-Protein Interaction Discovery via AlphaFold-Multimer*.

- LIA Threshold: ≥1610

- LIS Threshold: ≥0.073

We developed a small [script](../pipeline/lia_lis.py) in attaching average LIA/LIS scores from our AlphaPulldown workflow. 

**Workflow Image**
![Workflow](../images/workflow/variation_1.heic)

1. Extracting Features
```python
    excel_df = load_data(excel_path)
    excel_df['numerical'] = excel_df['jobs'].apply(extract_numeric) # Extract the numeric part from the 'jobs' column and create a new 'numerical' column
    excel_df_sorted = excel_df.sort_values(by='numerical') # Sort the DataFrame by the 'numerical' column
    excel_df_sorted = excel_df_sorted.drop(columns=['numerical']) # Drop the 'numerical' column as it is no longer needed
```

2. Filter for LIA/LIS
```python
def filter_features(data, thresholds):
    for feature, threshold in thresholds.items():
        data = data[data[feature] >= threshold]
    return data
...
    #filter data
    lia_threhsold = {
        "lis_score": 0.073,
        "average lia score": 1610
    }
    filter_df = filter_features(excel_df_sorted, lia_threhsold)

```

3. Finding Multipliers
```python
   #making a new column to grab the multiplier 
    excel_df_sorted['multiplier'] = excel_df_sorted['jobs'].apply(
    lambda x: -1 if x in filter_df['jobs'].values else 1)
    excel_df_sorted['iptm multiplied']=excel_df_sorted['iptm']*excel_df_sorted['multiplier']
    excel_df_sorted['mpdockq multiplied']=excel_df_sorted['mpDockQ/pDockQ']*excel_df_sorted['multiplier']
```

4. Multiplying constant and mapping metric to amino acid
```python
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
```

5. Determining DNFs

### Dependencies to download:
```bash
pip install pandas numpy absl-py biopython openpyxl
```

You can also [click here](LIA_LIS/variation_1_lialis.yml) to download the conda environment.

### How to run: 
This is the [code](LIA_LIS/variation1_lia.py) to run. 
```bash
python3 variation1_lia.py -uniprot=uniprot_id -file=excel_output_name -excel=/path/to/AlphaPulldown/outputs
```

Here is an example of the [output](../pipeline/example/flt3_variation1_lia_lis.xlsx)

---

## References
[1] Basu S, Wallner B. *DockQ: A Quality Measure for Protein-Protein Docking Models.* **PLOS ONE**, 11, 8, (2016). [Paper Link](https://doi.org/10.1371/journal.pone.0161879)

[2] Bertoline L., et al. *Before and after AlphaFold2: An Overview of Protein Structure Prediction*. **Frontiers in Bioinformatics**, 3, (2023). [Paper Link](https://doi.org/10.3389/fbinf.2023.1120370)

[3] Bryant P., Pozzati G., Zhu W. et al. *Predicting the structure of large protein complexes using AlphaFold and Monte Carlo tree search*. **Nat Commun**, 13, 6028 (2022). [Paper Link](https://doi.org/10.1038/s41467-022-33729-4)
  
[4] Kim AR, Hu Y, Comjean A, Rodiger J, Mohr SE, Perrimon N. *Enhanced Protein-Protein Interaction Discovery via AlphaFold-Multimer*. **BioRxiv** (2024). [Paper Link](https://www.biorxiv.org/content/10.1101/2024.02.19.580970v1)

[5] Xu J., Zhang Y. *How significant is a protein structure similarity with TM-score = 0.5?.* **Bioinformatics**, 26, 7, (2010). [Paper Link](https://doi.org/10.1093/bioinformatics/btq066)
