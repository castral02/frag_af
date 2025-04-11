# Variation 3
Variation 3 is a XGBoost Trained Model which classifies different fragments to No Pass, Low, Medium, or High Priority to discover dominant negative fragment. The data is based off of this [paper](https://www.cell.com/cell-systems/pdfExtended/S2405-4712(21)00157-5). 

We developed this depending on the previous models and what features are most prominent in similarity between each classification. We noticed that the most consistently accurate model is always filtered with mpDockQ/pDockQ score. From there, we grabbed features of ipTM and other interface structural features to help with the classifcation. 

## Dependencies to Download
To download the dependencies, [click here](xgboost.yml)

## The model itself:
### Overview
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
