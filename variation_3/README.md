# Variation 3: XGBoost Trained Model
Using the knowledge on Variation 2, we decided to create a trained model. Variation 2 at high accuracy was able to determine what is **not** a DNF but had a difficulty on determining what **is** a DNF. To increase precision, we developed a *XGBoost Trained Model.* XGBoost, Extreme Gradient Boosting, is a gradient boosted decision tree which is best used for regression, classification, and ranking problems (2). In our case, we aim to classify and rank the different fragments into **high, medium, and low** priority fragments. 

Before training the model, we had to clean our data (1). Our data is based of the Cell Systems Methods paper on Protein Tiling-- *Peptide-tiling screens of cancer drivers reveal oncogenic protein domains and associated peptide inhibitors* (3). However, due to the way we tiled the fragments, many of the DNFs are small interacting areas within the fragment itself making it difficult to categorize the fragment. In other words, the signal to noise ratio discrimination for this model to classify more difficult.

Classification of fragment is based on the percentage of similar amino acids to the original DNF fragment:

High: +95%

Medium: 75-95%

Low: 50-75%

With this in mind, our dataset was cleaned in two ways:
1. We took out proteins that have small regions of DNFs from the original dataset out
2. Fragments that score less than an average of 60% in Variation 1 and Variation 2

By cleaning the dataset, we hypothesized that the model will have a stronger predictive power in terms of accuracy and precision. To look at the trained dataset, [click here](library_dnf.csv)

Based on performances from Variations 1 and 2, mpDockQ, alone, had the best average accuracy; thus, we decided to use it as the threshold to lower the noise. 

![Accuracy Violin Plots]()

```python
# Filter features based on thresholds
def filter_features(data, thresholds):
    for feature, threshold in thresholds.items():
        data = data[data[feature] >= threshold]
    return data
...
thresholds = {'mpDockQ/pDockQ': 0.175
    }
filtered_data = filter_features(data, thresholds)
```

## Dependencies to Download
To download the dependencies, [click here](xgboost.yml)

## The model itself:
### Model 
This project implements an XGBoost-based machine learning model for classifying protein fragments based on their structural and interaction properties. The model processes protein fragment data with features such as polarity measurements, contact pairs, scoring metrics (sc), pi-scores, and AlphaFold confidence metrics (iptm, mpDockQ/pDockQ) to predict priority classes (Not Pass, Low, Medium, High).

### Key Features
- **Data Preprocessing**: Handles missing values, normalizes features per protein group, and applies feature filtering based on configurable thresholds
- **Class Imbalance Handling**: Implements SMOTE (Synthetic Minority Over-sampling Technique) to address class imbalance issues
- **Hyperparameter Optimization**: Uses GridSearchCV to find optimal XGBoost parameters
- **Performance Evaluation**: Provides comprehensive metrics including confusion matrices, accuracy, precision, recall, and F1 scores
- **Feature Selection**: Allows configurable feature selection for model training

### Model Performance
To look at the [log](training_info.log)
- **Best Hyperparameters**: 
  - colsample_bytree: 1.0
  - learning_rate: 0.1
  - max_depth: 5
  - n_estimators: 100
  - subsample: 1.0
- **Class Distribution**:
  - Original: Not Pass (0): 130, Low (1): 25, Medium (2): 8, High (3): 15
  - After SMOTE: Not Pass (0): 130, Low (1): 25, Medium (2): 130, High (3): 15
- **Accuracy**: 73%
- **Macro F1 Score**: 0.49
- **Per-Class Performance**:
  - Not Pass (0): Precision: 0.73, Recall: 0.62, F1: 0.67
  - Low (1): Precision: 0.00, Recall: 0.00, F1: 0.00
  - Medium (2): Precision: 0.84, Recall: 0.96, F1: 0.90
  - High (3): Precision: 0.50, Recall: 0.33, F1: 0.40

### Technical Details
- **Algorithm**: XGBoost with multiclass classification capability
- **Input Features**:
  - Polar: number of Polar Residues
  - contact_pairs: number of atomic contacts between residues 
  - sc: geometric shape commplementarity of protein-protein interface
  - [pi_score](https://www.nature.com/articles/s41467-021-23692-x): assessment of protein interaction 
  - iptm: AlphaFold confidence metric
  - mpDockQ/pDockQ: Docking quality scores
- **Output**: Classification into four categories (0: not_pass, 1: Low, 2: Medium, 3: High)
- **Performance Optimization**: Cross-validation with F1-macro scoring

### Notes on Model Limitations
- The model shows strong performance for Medium (class 2) predictions but struggles with Low (class 1) predictions
- Class imbalance remains a challenge despite SMOTE application
- Further optimization may be needed to improve performance on minority classes
    - It is not over representing high priority fragments; however, the model does know when the fragment is not a dominant negative fragment

### Usage
1. Configure threshold parameters for feature filtering
2. Select relevant features for model training
3. Run the script with the path to your CSV data file
4. The model will be trained, evaluated, and saved to the specified path

To train [model](model_creation.py), and here is the [training data set](library_dnf.csv)

To grab a trained [model](model_xgboost_moodel.pkl)

### Development Notes
The model normalizes features within each protein group and implements adaptive SMOTE parameters based on the smallest class size to ensure robust sampling even with limited data.

## How to run model:

The code to create your own predictions are based [here](prediction_dnf.py)

The way you run the code:
``` bash
python3 prediction_dnf.py -uniprot=uniprot_id -file=excel_output_name -excel=/path/to/AlphaPulldown/outputs
```

An example of the [output](../pipeline/example/flt3_version3_output.xlsx)

---
# References
[1] Ding, B., & Koutris, P., *Model selection for machine learning: The best choice is not always the most accurate.* **IEEE Data Engineering Bulletin,** 44(1), 24–33, (2021) [Paper Link](http://sites.computer.org/debull/A21mar/p24.pdf)

[2] NVIDIA, *XGBoost,* **NVIDIA,** (2025), [Link](https://www.nvidia.com/en-us/glossary/xgboost/)

[3] Ford K. et al., *Peptide-tiling screens of cancer drivers reveal oncogenic protein domains and associated peptide inhibitors,* 12(7), 716-732, (2021) [Paper Link](https://doi.org/10.1016/j.cels.2021.05.002)
