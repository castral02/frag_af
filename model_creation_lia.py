import pandas as pd
import numpy as np
import logging
from imblearn.over_sampling import SMOTE
from sklearn.cluster import Birch
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, GridSearchCV, cross_val_score
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report, precision_score, recall_score, f1_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
import joblib
import random
import os
import pickle

# Configure logging
logging.basicConfig(
    filename='/Users/castroverdeac/Desktop/codes_for_AF/actual_frag_af/model_variations/variation_3/RF/log/model_lia.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Load data from file
def load_data(filepath):
    return pd.read_excel(filepath)

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

# Perform BIRCH clustering
def birch_clustering(data):
    birch = Birch(branching_factor=50, threshold=0.3, n_clusters=4)
    labels = birch.fit_predict(data)
    return labels, birch

# Calculate BIRCH accuracy
def calculate_birch_accuracy(labels, ground_truth):
    return accuracy_score(ground_truth, labels)

# Apply SMOTE for oversampling
def oversample_with_smote(X, y):
    # Ensure there are enough samples to apply SMOTE
    if len(X) > 1:
        smote = SMOTE(sampling_strategy='auto', random_state=42, k_neighbors=1)
        return smote.fit_resample(X, y)
    else:
        # If not enough samples, return the original data without resampling
        logging.warning("Not enough samples to apply SMOTE. Returning original data.")
        return X, y


# Setup pipeline for scaling, clustering, and classification
# Log model performance metrics
def log_performance(metrics, params):
    logging.info(f"Model Parameters: {params}")
    logging.info(f"Confusion Matrix:\n{metrics['confusion_matrix']}")
    logging.info(f"Accuracy: {metrics['accuracy']:.2f}")
    logging.info(f"Precision: {metrics['precision']:.2f}")
    logging.info(f"Recall: {metrics['recall']:.2f}")
    logging.info(f"F1 Score: {metrics['f1']:.2f}")
    logging.info(f"Cross-validation scores: {metrics['cv_scores']}")
    logging.info(f"Classification Report:\n{metrics['classification_report']}")


# Iterative pipeline training
def iterative_birch_random_forest(X, y_true, rf_params, max_iterations=10, tolerance=0.01):
    best_combined_accuracy = 0
    best_rf_model = None
    updated_birch_params = {'threshold': 0.5}  # Example initial threshold for BIRCH

    rf_params = {
        'n_estimators': [100, 200],
        'max_depth': [10, 20, None],
        'min_samples_split': [2, 5],
        'min_samples_leaf': [1, 2],
        'max_features': [None, 'sqrt']
    }

    for iteration in range(max_iterations):
        print(f"Iteration {iteration + 1}")
        
        # Step 1: BIRCH Clustering
        birch = Birch(**updated_birch_params)
        birch.fit(X)
        y_pred_birch = birch.labels_
        
        # Calculate BIRCH accuracy
        birch_accuracy = accuracy_score(y_true, y_pred_birch)
        print(f"BIRCH Accuracy: {birch_accuracy}")
        
        # Step 2: Train Random Forest with SMOTE
        X_train, X_test, y_train, y_test = train_test_split(X, y_true, test_size=0.2, random_state=42)
        
        # Apply SMOTE on training data
        smote = SMOTE(random_state=42, k_neighbors=1)
        X_train_resampled, y_train_resampled = smote.fit_resample(X_train, y_train)
        print(f"Resampled training data: {len(X_train_resampled)} samples")

        rf = RandomForestClassifier(random_state=42)
        grid_search = GridSearchCV(
            estimator=rf,
            param_grid=rf_params,
            scoring='accuracy',
            cv=3
        )
        grid_search.fit(X_train_resampled, y_train_resampled)
        
        best_rf = grid_search.best_estimator_
        rf_accuracy = grid_search.best_score_
        print(f"Random Forest Accuracy: {rf_accuracy}")
        
        # Step 3: Calculate Combined Score
        combined_accuracy = 0.5 * birch_accuracy + 0.5 * rf_accuracy
        print(f"Combined Accuracy: {combined_accuracy}")
        
        # Check for improvement
        if combined_accuracy - best_combined_accuracy < tolerance:
            print("Stopping criteria met: Minimal improvement in accuracy.")
            break
        
        best_combined_accuracy = combined_accuracy
        best_rf_model = best_rf
        
        # Step 4: Update BIRCH Parameters
        if birch_accuracy < rf_accuracy:
            updated_birch_params['threshold'] = max(0.1, updated_birch_params['threshold'] - 0.05)
        else:
            updated_birch_params['threshold'] = min(1.0, updated_birch_params['threshold'] + 0.05)
        

        print(f"Updated BIRCH Parameters: {updated_birch_params}\n")
    
    # Final Evaluation
    print(f"Best Combined Accuracy Achieved: {best_combined_accuracy}")
    y_test_pred = best_rf_model.predict(X_test)
    
    # Compute the metrics
    accuracy = accuracy_score(y_test, y_test_pred)
    precision = precision_score(y_test, y_test_pred, average='weighted')
    recall = recall_score(y_test, y_test_pred, average='weighted')
    f1 = f1_score(y_test, y_test_pred, average='weighted')
    confusion = confusion_matrix(y_test, y_test_pred)
    cv_scores = cross_val_score(best_rf_model, X, y_true, cv=3)
    class_report = classification_report(
        y_test, y_test_pred, target_names=['not_pass', 'Low', 'Medium', 'High'], labels=[0, 1, 2, 3]
    )    
    # Prepare the metrics dictionary
    metrics = {
        'accuracy': accuracy,
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'confusion_matrix': confusion,
        'cv_scores': cv_scores,
        'classification_report': class_report
    }
    
    # Log the performance
    log_performance(metrics, {'BIRCH': updated_birch_params, 'RandomForest': rf_params})
    
    return best_combined_accuracy, best_rf_model, updated_birch_params


# Main function to tie everything together
def main(filepath, thresholds, rf_iterations=10):
    # Load and preprocess data
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")

    data = load_data(filepath)
    data = preprocess_labels(data)  # Preprocess labels before filtering
    
    # Apply feature filtering
    filtered_data = filter_features(data, thresholds)
    y_true = filtered_data["known_label"]
    
    # Scale features for clustering
    cluster_features = ['Polar', 'Hydrophobhic', 'sc', 'pi_score', 'iptm', 'average pae score']
    X_cluster = filtered_data[cluster_features]
    scaler = StandardScaler()
    X_cluster_scaled = scaler.fit_transform(X_cluster)

    # Perform BIRCH clustering
    birch_labels, _ = birch_clustering(X_cluster_scaled)
    birch_accuracy = calculate_birch_accuracy(birch_labels, y_true)
    logging.info(f"BIRCH Clustering Accuracy: {birch_accuracy:.2f}")
    
    # Train Random Forest with iterative refinement
    rf_accuracy, best_model, updated_birch_params = iterative_birch_random_forest(
        X_cluster_scaled, y_true, rf_iterations
    )
    logging.info(f"Best Random Forest Accuracy: {rf_accuracy:.2f}")
    
    # Evaluate model performance
    y_pred = best_model.predict(X_cluster_scaled)
    accuracy = accuracy_score(y_true, y_pred)
    cm = confusion_matrix(y_true, y_pred)
    logging.info(f"Confusion Matrix:\n{cm}")
    logging.info(f"Final Accuracy: {accuracy:.2f}")
    
    return best_model

# Example usage
if __name__ == "__main__":
    thresholds = {
        "lis_score": 0.073,
        "average lia score": 1610
    }
    filepath = "/Users/castroverdeac/Desktop/codes_for_AF/actual_frag_af/data/data_set_correct.xlsx"
    best_model=main(filepath, thresholds)
    # Example: Assuming `model` is your trained model
    filename = '/Users/castroverdeac/Desktop/codes_for_AF/actual_frag_af/model_variations/variation_3/RF/models/model_lia.pkl'

# Save the model to a file
    with open(filename, 'wb') as file:
        pickle.dump(best_model, file)
