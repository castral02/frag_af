# Frag-AF
A trained clustering model in discovering small inhibitory peptide binders, also known as dominant negative fragments, for a high throughput method in understanding protein function. 

## Abstract
I can help you develop this abstract further. Here's an expanded version that includes your model testing on metabolic proteins:

# Abstract

Peptide tiling is a method in which small segments of a protein interactor are screened to identify fragments capable of exerting a regulatory function, such as inhibition. One of the challenges of this approach is that—due to the vast size of the protein sequence space—screening large libraries of fragments for potential interactors is often prohibitively costly in terms of time and resources. New protein structure prediction technologies, such as Google DeepMind's AlphaFold, have the potential to decrease experimental hurdles by allowing researchers to focus on highly confident interacting areas.

In this study, we developed an XGBoost-based classification model to identify dominant negative protein fragments with regulatory potential. The model was trained on features including number of polar residues at the interface, contact pairs, scoring coefficients, and AlphaFold confidence metrics. We applied Synthetic Minority Over-sampling Technique (SMOTE) to address class imbalance challenges and optimized model parameters through grid search cross-validation. 

The resulting classifier achieved 73% overall accuracy with particularly strong performance in identifying medium-priority fragments (F1-score: 0.90). We validated our approach by testing the model on 50 novel metabolic proteins, demonstrating its potential to accelerate the discovery of bioactive peptide fragments while reducing experimental costs. This computational pre-screening tool represents a significant advancement in rational protein fragment design, offering researchers a method to prioritize candidates for experimental validation in therapeutic development pipelines.


## What are the different variations?
Mulitple Variations were tested and different metrics were tested in exploring these dominant negative fragments. To validate and test different metrics, we used previous [data](https://www.cell.com/cell-systems/pdfExtended/S2405-4712(21)00157-5). 

When doing the screening/peptide tiling, we utilized a looping iteration of AlphaFold Multimer called [AlphaPulldown](https://github.com/KosinskiLab/AlphaPulldown). We used version 1.0.4 off of Biowulf. For MSA creation, we used a ColabFold Search. 

[Click Here](sbatch_files_example/) to see AlphaPulldown pipeline.  

### Variation 1
In this variation, we correlated previously known metrics to previous data. We explored four metrics: Local Interaction Area, Local Interaction Score, mpDockQ, and ipTM. 

[Click Here](variation_1) to get the codes and information

### Variation 2
In this variation, we used acceptable thresholding metrics to filter fragments and then cluster the items based on different structural and AlphaFold metrics. There are several avenues we took in looking at this. We decided to go single thresholds of just ipTM and mpDockQ or dual thresholds were you look at LIA and LIS or mpDockQ and ipTM. We saw that there is some predictive power in the metrics; however, the precision is not that strong, therefore, we believed we can cluster the data to increase precision. 

[Click Here](variation_2) to get the codes and information. 

### Variation 3 
In this variation, we created a trained classification machine learning model based on previous data. According to previous variations and models, the best metric for quantification is mpDockQ alone with a threshold in terms of accuracy and percision for a multitiude of proteins. 

[Click Here](variation_3) to get the codes and information. 

## Declaration of generative AI usage
This project utilized OpenAI's ChatGPT to assist in generating Python code, documentation, or other textual content.

## Citation
Biowulf Acknowledgement: This work utilized the computational resources of the NIH HPC [Biowulf cluster](https://hpc.nih.gov).

AlphaPulldown: Dingquan Yu, Grzegorz Chojnowski, Maria Rosenthal, Jan Kosinski, AlphaPulldown—a python package for protein–protein interaction screens using AlphaFold-Multimer, Bioinformatics, Volume 39, Issue 1, January 2023, btac749, [paper link](https://doi.org/10.1093/bioinformatics/btac749)

LIA/LIS: Kim AR, Hu Y, Comjean A, Rodiger J, Mohr SE, Perrimon N. "Enhanced Protein-Protein Interaction Discovery via AlphaFold-Multimer" bioRxiv (2024) [paper link](https://www.biorxiv.org/content/10.1101/2024.02.19.580970v1)
