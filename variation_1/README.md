# Variation 1
As initial approaches, we correlated commonly used AlphaFold Metrics to dominant negative fragments. In this repository, we explored 3 metircs: Local Interaction Area and Local Interaction Score, mpDockQ, and ipTM. The workflow for each of these codes stay the same throughout.

1. AlphaPulldown Workflow
2. Extracting AlphaFold Metrics
3. Mapping Metrics onto individual amino acids
4. Normalizing Metrics
5. Subtracting experimentally known thresholds to the metrics
6. Inferring Dominant Negative Fragments (DNFs)

**Workflow Image**

![Workflow](../images/line_plots/workflow/frag_af_variation_1.heic)

## LIA/LIS
Local Interaction Area and Local Interaction Score is a newly developed metric from the Perrimon Lab at Harvard. The scoring system is derived off of the AlphaFold Metric Predicted Aligned Error to discover highly interactive areas when the strucutre is small and felxible.

We developed a small [script](../pipeline/lia_lis.py) in attaching average LIA/LIS scores from our AlphaPulldown workflow. 

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

---

## References
- Bertoline, Letícia M. F., et al. *Before and after AlphaFold2: An Overview of Protein Structure Prediction*. **Frontiers in Bioinformatics**, 3, (2023). [Paper Link](https://doi.org/10.3389/fbinf.2023.1120370)

- Bryant, P., Pozzati, G., Zhu, W. et al. *Predicting the structure of large protein complexes using AlphaFold and Monte Carlo tree search*. **Nat Commun**, 13, 6028 (2022). [Paper Link](https://doi.org/10.1038/s41467-022-33729-4)
- 
- Kim AR, Hu Y, Comjean A, Rodiger J, Mohr SE, Perrimon N. *Enhanced Protein-Protein Interaction Discovery via AlphaFold-Multimer*. **bioRxiv** (2024). [Paper Link](https://www.biorxiv.org/content/10.1101/2024.02.19.580970v1)
