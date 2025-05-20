# Variation 2: Unsupervised BIRCH Clustering
We evolved variation 1 into an unsupervised clustering framework. We generated 60 amino acid residue fragments using a sliding window of 10 residues, resulting in a 50-residue overlap between adjacent fragments. We hypothesized that the 50 amino acid overlap would cause fragments to have similar interface features and AlphaFold Metrics (i.e. mpDockQ and ipTM) allowing us to cluster fragments due to their comparable binding modes. To capture these relationships, we applied **BIRCH (Balanced Iterative Reducing and Clustering using Hierarchies)**. BIRCH is a hierarchical clustering algorithm that is well-suited for large datasets, and for our case, incremental learning (1).

We used the BIRCH system due to its ability to represent fragment similarity in a hierarchical branching structure, allowing an intuitive classification into **high, medium, and low** priority clusters based on interaction potentials. The system reflects the underlying biological assumption that protein fragments with similar sequences will have similar interface characteristics and therefore exhibit similar binding modes (2).

We implemented two strategies: 

1. Single Threshold: mpDockQ and ipTM allowing us to assess the individual contribution of each metric separately
```python
mpdockq_threhsold = {
        "mpDockQ/pDockQ": 0.175}
iptm_threhshold = {
       "iptm": 0.5}
```

2. Dual Threshold: LIA/LIS together or mpDockQ with ipTM allowing a holistic view of protein-protein interaction quality 
```python
lia_threhsold = {
        "lis_score": 0.073,
        "average lia score": 1610}
```

```python
mpdock_threhshold = {
      "mpDockQ/pDockQ": 0.175,
       "iptm": 0.5}
```

**Workflow**
![Workflow Image](../images/workflow/variation_2.heic)

1. Extracting data
```python
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
```

2. Filtering Fragments
```python
def filter_features(data, thresholds):
    for feature, threshold in thresholds.items():
        data = data[data[feature] >= threshold]
    return data
...
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
```

3. Normalizing data
```python
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
```

4. Clustering data
```python
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
```

5. Mapping data to individual amino acids
```python
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
```
6. Determining DNFs

## Single Threshold
### Dependencies to download
``` bash
pip install absl-py matplotlib numpy pandas requests scikit-learn seaborn openpyxl
```
Or you can [click here](single_threshold/variation_2_single_threshold.yml )to download the conda environment 

### How to run:
This is the [code](single_threshold/variation2_single_thresholds.py)

```bash
python3 variation2_single_thresholds.py -uniprot=uniprot_id -file=excel_output_name -excel=/path/to/AlphaPulldown/outputs
```

Here is an [example](../pipeline/example/flt3_single_threshold.xlsx) of an output.

## Dual Threshold
### Dependencies to download
``` bash
pip install absl-py matplotlib numpy pandas requests scikit-learn seaborn openpyxl
```

Or you can [click here](dual_threshold/variation_2_dual_threshold.yml )to download the conda environment 

### How to run:
This is the [code](dual_threshold/variation_2.py)

```bash
python3 variation_2.py -uniprot=uniprot_id -file=excel_output_name -excel=/path/to/AlphaPulldown/outputs
```

Here is an [example](../pipeline/example/flt3_dual_threshold_v2.xlsx) of an output.

----

## References
[1] Zhang T., Ramakrishnan R., & Livny M., *BIRCH: an efficient data clustering method for very large databases,* **ACM Sigmod Record,** 25(2), (1996) [Paper Link](https://doi.org/10.1145/235968.233324)

[2] Madaoui H., & Guerois R., *Coevolution at protein complex interfaces can be detected by the complementarity trace with important impact for predictive docking,* **Proc. Natl. Acad. Sci. U.S.A.,** 105(22), (2008), [Paper Link](https://doi.org/10.1073/pnas.0707032105)
