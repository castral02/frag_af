# Variation 2
We evolved variation 1 into an unsupervised clustering system assuming that similar fragments will have similar features and metrics. We hypothesized that system would increase the accuracy and predictive power of inferring these DNFs. We developed two pipelines: single threshold versus dual threshold. Single threshold is exploring ipTM and mpDockQ independently; while dual is exploring ipTM and mpDockQ together. We extended this dual threshold theme using LIA/LIS which is the preferred method when determing protein-protein interaction. 

**Workflow**
1. Extracting data
2. Filtering Fragments
3. Normalizing data
4. Clustering data
5. Mapping data to individual amino acids
6. Determining DNFs

![Workflow Image]()

## Single Threshold
In the single threshold, we only did ipTM and mpDockQ separately. 
```bash
lia_threhsold = {
        "mpDockQ/pDockQ": 0.175}
mpdock_threhshold = {
       "iptm": 0.5}
```

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
For Dual Threshold, there were two types of filtering systems we did:
1. LIA/LIS
```bash
lia_threhsold = {
        "lis_score": 0.073,
        "average lia score": 1610}
```

2. mpDockQ and ipTM
```bash
mpdock_threhshold = {
      "mpDockQ/pDockQ": 0.175,
       "iptm": 0.5}
```

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
