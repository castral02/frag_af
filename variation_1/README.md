# Variation 1
In this variation, we correlated previously known metrics to previous data. We explored four metrics: Local Interaction Area, Local Interaction Score, mpDockQ, and ipTM. 

## LIA/LIS
[LIA and LIS](https://github.com/flyark/AFM-LIS) is a newer metric in exploring protein-protein interactions when there is not enough predictive power. In our case, we continued using their average threshold metrics. 

In order to grab LIA/LIS, we used their code to grab the materials from our AlphaPulldown and attach the LIA/LIS to the AlphaPulldown CSV. 

Use this [code](../pipeline/lia_lis.py) to grab the LIA and LIS and attach it to the csv. 

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

## mpDockQ/ipTM
In this code, we decided to use two different metrics that are popular to explore [ipTM](https://academic.oup.com/bioinformatics/article/26/7/889/213219?login=true) and [mpDockQ](https://www.nature.com/articles/s41467-022-33729-4).

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
