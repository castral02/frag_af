# Variation 2
In this variation, we used acceptable thresholding metrics to filter fragments and then cluster the items based on different structural and AlphaFold metrics. There are several avenues we took in looking at this. We decided to go single thresholds of just ipTM and mpDockQ or dual thresholds were you look at LIA and LIS or mpDockQ and ipTM.
We saw that there is some predictive power in the metrics from version 1; however, the precision is not that strong, therefore, we believed we can cluster the data to increase precision.

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
