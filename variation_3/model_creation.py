import pandas as pd
import numpy as np
import logging
from imblearn.over_sampling import SMOTE
from sklearn.model_selection import train_test_split, GridSearchCV, cross_val_score
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report, precision_score, recall_score, f1_score
import os
import pickle
import xgboost as xgb
from sklearn.preprocessing import MinMaxScaler
from collections import Counter


logging.basicConfig(
    filename='/Users/castroverdeac/Desktop/codes_for_AF/actual_frag_af/model_variations/variation_3/XGBoost/log/model_xgboost_mpdock.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Load data from file
def load_data(filepath):
    return pd.read_csv(filepath)

# Filter features based on thresholds
def filter_features(data, thresholds):
    for feature, threshold in thresholds.items():
        data = data[data[feature] >= threshold]
    return data

# Preprocess priority labels
def preprocess_labels(data):
    label_mapping = {'not_pass': 0, 'Low': 1, 'Medium': 2, 'High': 3}
    data['known_label'] = data['known_label'].map(label_mapping).fillna(0).astype(int)
    return data

# Normalize Features for Each Protein and Dropping the Nones
def normalize_drop(data):
    # Remove stray spaces from column names
    data.columns = data.columns.str.strip()

    # Drop rows with NaNs in numeric columns
    numeric_cols = ['Polar',
                         'contact_pairs','sc',
                         'pi_score',
                         'mpDockQ/pDockQ']
    data = data.dropna(subset=numeric_cols)

    # Exclude label column from normalization
    numeric_cols = [col for col in numeric_cols if col != 'known_label']

    # Normalize per tiled_protein
    data[numeric_cols] = data.groupby('tiled protein')[numeric_cols].transform(
        lambda x: MinMaxScaler().fit_transform(x.values.reshape(-1, 1)).flatten()
    )
    return data

# Apply SMOTE for oversampling
def oversample_with_smote(X, y, strategy='minority'):
    """
    Apply SMOTE to oversample underrepresented classes.
    
    Args:
        X (array-like): Feature matrix.
        y (array-like): Labels.
        strategy (str or dict): SMOTE sampling strategy. Default 'auto' = balance all classes.
    
    Returns:
        X_resampled, y_resampled
    """
    if len(X) <= 1:
        logging.warning("Not enough samples to apply SMOTE. Returning original data.")
        return X, y

    # Count original class distribution
    class_counts = Counter(y)
    logging.info(f"Original class distribution: {class_counts}")

    try:
        # Adjust k_neighbors based on smallest class
        min_class_size = min(class_counts.values())
        k_neighbors = min(5, min_class_size - 1) if min_class_size > 1 else 1

        smote = SMOTE(sampling_strategy=strategy, random_state=42, k_neighbors=k_neighbors)
        X_resampled, y_resampled = smote.fit_resample(X, y)

        new_counts = Counter(y_resampled)
        logging.info(f"Resampled class distribution: {new_counts}")

        return X_resampled, y_resampled

    except ValueError as e:
        logging.error(f"SMOTE failed: {e}")
        return X, y

def log_performance(model, X_test, y_test, params, cv_scores=None):
    # Make predictions
    y_pred = model.predict(X_test)
    
    # Calculate metrics
    confusion = confusion_matrix(y_test, y_pred)
    accuracy = accuracy_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred, average='macro', zero_division=0)
    recall = recall_score(y_test, y_pred, average='macro', zero_division=0)
    f1 = f1_score(y_test, y_pred, average='macro', zero_division=0)
    report = classification_report(y_test, y_pred)
    
    # Create metrics dictionary
    metrics = {
        'confusion_matrix': confusion,
        'accuracy': accuracy,
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'classification_report': report,
        'cv_scores': cv_scores if cv_scores is not None and len(cv_scores) > 0 else "N/A"
    }
    
    # Log the performance
    logging.info(f"Model Parameters: {params}")
    logging.info(f"Confusion Matrix:\n{metrics['confusion_matrix']}")
    logging.info(f"Accuracy: {metrics['accuracy']:.2f}")
    logging.info(f"Precision: {metrics['precision']:.2f}")
    logging.info(f"Recall: {metrics['recall']:.2f}")
    logging.info(f"F1 Score: {metrics['f1']:.2f}")
    logging.info(f"Cross-validation scores: {metrics['cv_scores']}")
    logging.info(f"Classification Report:\n{metrics['classification_report']}")

def select_relevant_features(data, selected_features):
    """
    Filters the data to only include the selected features.
    """
    return data[selected_features]

def main(filepath, thresholds, selected_features, rf_iterations=10):
    # Load data
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    
    data = load_data(filepath)
    data = preprocess_labels(data)  # Preprocess labels before filtering

    # Filter the data based on the feature thresholds
    filtered_data = filter_features(data, thresholds)

    # Normalize data
    filtered_data=normalize_drop(filtered_data)

    # Select only the relevant features for training
    filtered_data = select_relevant_features(filtered_data, selected_features)
    
    X = filtered_data.drop(columns=['known_label'])
    y = filtered_data['known_label']
    
    # Apply SMOTE for oversampling
    X_train_resampled, y_train_resampled = oversample_with_smote(X, y)

    # Now split the resampled training data into new training and test sets
    X_train, X_test, y_train, y_test = train_test_split(X_train_resampled, y_train_resampled, test_size=0.2, random_state=42)

    # Initialize XGBoost model
    xgboost_model = xgb.XGBClassifier(
        random_state=42,
        eval_metric='mlogloss',
        use_label_encoder=False
    )

    # Train the model
    # Define parameter grid for GridSearchCV
    param_grid = {
        'max_depth': [3, 5, 7],
        'learning_rate': [0.01, 0.1, 0.2],
        'n_estimators': [100, 200, 300],
        'subsample': [0.8, 1.0],
        'colsample_bytree': [0.8, 1.0]
    }

    # Grid Search for tuning XGBoost hyperparameters
    grid_search = GridSearchCV(estimator=xgboost_model, param_grid=param_grid, cv=5, scoring='f1_macro', n_jobs=-1)
    
    grid_search.fit(X_train, y_train)

    # Best hyperparameters
    logging.info(f"Best Parameters: {grid_search.best_params_}")

    # Evaluate the model with best hyperparameters
    best_model = grid_search.best_estimator_
    y_pred_best = best_model.predict(X_test)

    # Log evaluation metrics
    log_performance(
        model=best_model,
        X_test=X_test,
        y_test=y_test,
        params=grid_search.best_params_,
        cv_scores=grid_search.cv_results_['mean_test_score']
    )

    # Save the best model to a file
    filename = '/Users/castroverdeac/Desktop/codes_for_AF/actual_frag_af/model_variations/variation_3/XGBoost/model/model_xgboost_model.pkl'
    with open(filename, 'wb') as file:
        pickle.dump(best_model, file)

if __name__ == "__main__":
    thresholds = {'mpDockQ/pDockQ': 0.175
    }

    selected_features = ['Polar',
                         'contact_pairs','sc',
                         'pi_score',
                         'iptm','known_label']  

    filepath = "/Users/castroverdeac/Desktop/codes_for_AF/actual_frag_af/model_variations/variation_3/XGBoost/library_dnf.csv"
    main(filepath, thresholds, selected_features)
