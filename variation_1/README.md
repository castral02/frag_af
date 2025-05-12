# Variation 1

As initial approaches, we correlated commonly used AlphaFold Metrics to dominant negative fragments. In this repository, we explored 3 metircs: Local Interaction Area and Local Interaction Score, mpDockQ, and ipTM. The workflow for each of these codes stay the same throughout.

1. AlphaPulldown Workflow
2. Extracting AlphaFold Metrics
3. Mapping Metrics onto individual amino acids
4. Normalizing Metrics
5. Subtracting experimentally known thresholds to the metrics
6. Inferring Dominant Negative Fragments (DNFs)

**Workflow Image**

![Workflow](../images/workflow/frag_af_variation_1.heic)

## LIA/LIS

Local Interaction Area (LIA) and Local Interaction Score (LIS) is a newly developed metric from the Perrimon Lab at Harvard. The scoring system is derived off of the AlphaFold Metric Predicted Aligned Error (PAE )to discover highly interactive areas when the strucutre is small and felxible (4). 

We used the average LIA/LIS scores of each fragment and the threshold's from the original paper-- *Enhanced Protein-Protein Interaction Discovery via AlphaFold-Multimer*.

- LIA Threshold: ≥1610

- LIS Threshold: ≥0.073

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

The other half of this variation is exploring two popular metrics: interface predicted templating modelling ([ipTM](https://academic.oup.com/bioinformatics/article/26/7/889/213219?login=true)) and [mpDockQ](https://www.nature.com/articles/s41467-022-33729-4).

Thresholds were determined from previous research:

- mpDockQ: ≥ 0.173 (1, 3)
- ipTM: ≥ 0.5 (2, 5)

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
[1] Basu S, Wallner B. *DockQ: A Quality Measure for Protein-Protein Docking Models.* **PLOS ONE**, 11, 8, (2016). [Paper Link](https://doi.org/10.1371/journal.pone.0161879)

[2] Bertoline L., et al. *Before and after AlphaFold2: An Overview of Protein Structure Prediction*. **Frontiers in Bioinformatics**, 3, (2023). [Paper Link](https://doi.org/10.3389/fbinf.2023.1120370)

[3] Bryant P., Pozzati G., Zhu W. et al. *Predicting the structure of large protein complexes using AlphaFold and Monte Carlo tree search*. **Nat Commun**, 13, 6028 (2022). [Paper Link](https://doi.org/10.1038/s41467-022-33729-4)
  
[4] Kim AR, Hu Y, Comjean A, Rodiger J, Mohr SE, Perrimon N. *Enhanced Protein-Protein Interaction Discovery via AlphaFold-Multimer*. **BioRxiv** (2024). [Paper Link](https://www.biorxiv.org/content/10.1101/2024.02.19.580970v1)

[5] Xu J., Zhang Y. *How significant is a protein structure similarity with TM-score = 0.5?.* **Bioinformatics**, 26, 7, (2010). [Paper Link](https://doi.org/10.1093/bioinformatics/btq066)
