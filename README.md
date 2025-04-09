# Frag-AF
A trained clustering model in discovering small inhibitory peptide binders, also known as dominant negative fragments, for a high throughput method in understanding protein function. 

## Abstract
Peptide tiling is a method in which small segments of a protein interactor are screened to identify fragments capable of exterting a regulatory function, for instance, inhibition. One of the challenges of this approach is that-- due to the vast size fo the protein sequence space-- screening large libraries of fragments for potentital interactors are often prohibitively costly in terms of time and resources. New protein structure prediction technologies, such as Google DeepMind's Alphafold, have the potential to decrease experimental hurdles by allowing researchers to focus on highly confident interacting areas. 

In this study, we developed a trained classifing model in discovering dominant negative fragments, and tested the model with metabolic proteins. 

## What are the different variations?
Mulitple Variations were tested and different metrics were tested in exploring these dominant negative fragments. To validate and test different metrics, we used previous [data](https://www.cell.com/cell-systems/pdfExtended/S2405-4712(21)00157-5). 

When doing the screening/peptide tiling, we utilized a looping iteration of AlphaFold Multimer called [AlphaPulldown](https://github.com/KosinskiLab/AlphaPulldown). We used version 1.0.4 off of Biowulf. For MSA creation, we used a ColabFold Search. 

[Click Here](sbatch_files_examples) to see AlphaPulldown pipeline.  

### Variation 1
In this variation, we correlated previously known metrics to previous data. We explored four metrics: Local Interaction Area, Local Interaction Score, mpDockQ, and ipTM. 

#### LIA/LIS
[LIA and LIS](https://github.com/flyark/AFM-LIS) is a newer metric in exploring protein-protein interactions when there is not enough predictive power. In our case, we continued using their average threshold metrics. 

In order to grab LIA/LIS, we used their code to grab the materials from our AlphaPulldown and attach the LIA/LIS to the AlphaPulldown CSV. 

Use this [code](pipeline/lia_lis.py) to grab the LIA and LIS and attach it to the csv. 

###### Dependencies to download:
```bash
pip install pandas numpy absl-py biopython
```

###### How to run: 
```bash
python3 lia_lis.py -output_dir /path/to/AlphaPulldown/materials
```
##### How to run LIA/LIS


## How to run the latest variation?

## Declaration of generative AI usage
This project utilized OpenAI's ChatGPT to assist in generating Python code, documentation, or other textual content.

## Citation
Biowulf Acknowledgement: This work utilized the computational resources of the NIH HPC [Biowulf cluster](https://hpc.nih.gov).
AlphaPulldown: Dingquan Yu, Grzegorz Chojnowski, Maria Rosenthal, Jan Kosinski, AlphaPulldown—a python package for protein–protein interaction screens using AlphaFold-Multimer, Bioinformatics, Volume 39, Issue 1, January 2023, btac749, [paper link](https://doi.org/10.1093/bioinformatics/btac749)
LIA/LIS: Kim AR, Hu Y, Comjean A, Rodiger J, Mohr SE, Perrimon N. "Enhanced Protein-Protein Interaction Discovery via AlphaFold-Multimer" bioRxiv (2024) [paper link](https://www.biorxiv.org/content/10.1101/2024.02.19.580970v1)
