# Imports

import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix, matthews_corrcoef
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from xgboost import XGBClassifier
from rnanorm import CPM  # Install with: pip install rnanorm

# ----------------------- 1. Preprocessing -----------------------

# This function takes a DataFrame, a label column name, and optional scaler and fit_scaler parameters.
# It returns the preprocessed features, encoded labels, scaler, and label encoder.
def preprocess_data(df, label_col='cases.disease_type', scaler=None, fit_scaler=True):
    y = df[label_col]
    X = df.drop(columns=[label_col])

    # Convert to CPM (Counts Per Million) to normalize library size
    X_normalized = CPM().fit_transform(X)

    X_log = np.log1p(X_normalized)  # Log transformation

    # Standardize the data
    if scaler is None:
        scaler = StandardScaler()
    if fit_scaler:
        X_scaled = scaler.fit_transform(X_log)
    else:
        X_scaled = scaler.transform(X_log)

    label_encoder = LabelEncoder()
    y_encoded = label_encoder.fit_transform(y)

    return X_scaled, y_encoded, scaler, label_encoder

# ----------------------- 2. Feature Selection -----------------------
# This function takes a DataFrame and a list of significant genes.
# It returns a list of significant genes that are present in the DataFrame.
def select_significant_genes(df, gene_list):
    all_genes = df.columns.drop('cases.disease_type')
    return [gene for gene in gene_list if gene in all_genes]

# ----------------------- 3. Training & Evaluation -----------------------
# This function takes a model name, model instance, training and testing data, label encoder, and optional feature names.
# It evaluates the model using cross-validation and computes various metrics.
def evaluate_model_sig(name, model, X_train, y_train, X_test, y_test, label_encoder, feature_names=None):
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    y_pred_cv = cross_val_predict(model, X_train, y_train, cv=skf)
    model.fit(X_train, y_train)
    y_pred_test = model.predict(X_test)

    def get_metrics(y_true, y_pred):
        return {
            'Accuracy': accuracy_score(y_true, y_pred),
            'Precision': precision_score(y_true, y_pred, average='weighted', zero_division=0),
            'Recall': recall_score(y_true, y_pred, average='weighted', zero_division=0),
            'F1': f1_score(y_true, y_pred, average='weighted', zero_division=0),
            'MCC': matthews_corrcoef(y_true, y_pred),
            'ConfusionMatrix': confusion_matrix(y_true, y_pred)
        }

    results = {
        'CV': get_metrics(y_train, y_pred_cv),
        'Test': get_metrics(y_test, y_pred_test)
    }
    # print the accuracy scores for each model
    print(f"\n{name} - Metrics:")
    print(f"CV Accuracy: {results['CV']['Accuracy']:.4f}")
    print(f"Test Accuracy: {results['Test']['Accuracy']:.4f}")
    

    cm = results['Test']['ConfusionMatrix']
    decoded_labels = label_encoder.inverse_transform(np.unique(y_test))
    plt.figure(figsize=(6, 4))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues',
                xticklabels=decoded_labels, yticklabels=decoded_labels)
    plt.title(f'{name} - Confusion Matrix (Test Set)')
    plt.xlabel('Predicted')
    plt.ylabel('Actual')
    plt.tight_layout()
    plt.show()


    top_features = []
    # Feature Importance Extraction
    if hasattr(model, 'feature_importances_') and feature_names is not None:
        importances = model.feature_importances_
        indices = np.argsort(importances)[::-1]  # Full descending sort
        top_features = np.array(feature_names)[indices].tolist()
        
    elif hasattr(model, 'coef_') and feature_names is not None:
        importances = np.abs(model.coef_).sum(axis=0)
        indices = np.argsort(importances)[::-1]  # Full descending sort
        top_features = np.array(feature_names)[indices].tolist()

    # After computing importances and indices:
    if feature_names is not None and importances is not None:
        # Get top 10
        top_n = 10
        top_indices = indices[:top_n]
        plt.figure(figsize=(8, 6))
        sns.barplot(
            x=importances[top_indices],
            y=np.array(feature_names)[top_indices],
            orient='h'
        )
        plt.title(f'{name} - Top 10 Feature Importances')
        plt.xlabel('Importance')
        plt.ylabel('Feature')
        plt.tight_layout()
        plt.show()


    return results, top_features 

# ----------------------- 4. Model Training -----------------------

# This function takes training and testing data, label encoder, gene set label, and feature names.
# It trains and evaluates multiple models (SVM, Logistic Regression, XGBoost) and returns their results and top genes.
def train_and_evaluate_models_sig(X_train, y_train, X_test, y_test, label_encoder, gene_set_label, feature_names):
    models = {
        'SVM': SVC(kernel='linear', probability=True),
        'LogisticRegression': LogisticRegression(max_iter=1000),
        'XGBoost': XGBClassifier(n_estimators=100, learning_rate=0.1, max_depth=5,
                                 random_state=42, eval_metric='mlogloss')
    }

    all_results = {}
    model_top_genes = {}
    for name, model in models.items():
        print(f"Training {name} with {gene_set_label}...")
        results, top_genes = evaluate_model_sig(name, model, X_train, y_train, X_test, y_test, label_encoder, feature_names)
        all_results[f"{name}_{gene_set_label}"] = results['Test']
        model_top_genes[f"{name}_{gene_set_label}"] = top_genes

    return all_results, model_top_genes

# ----------------------- 5. Plot Comparison -----------------------

# This function takes a dictionary of results for significant genes and plots a comparison of the metrics.
def plot_comparison(results_significant):
    metrics = ['Accuracy', 'Precision', 'Recall', 'F1', 'MCC']
    df = []

    for model_label, metrics_dict in results_significant.items():
        for metric in metrics:
            df.append({
                'Model': model_label,
                'Metric': metric,
                'Score': metrics_dict[metric]
            })

    df = pd.DataFrame(df)
    plt.figure(figsize=(8, 4))
    sns.barplot(data=df, x='Metric', y='Score', hue='Model')
    plt.title('Model Comparison Across Metrics')
    plt.xticks(rotation=45)
    plt.ylim(0.80, 1.00)
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.show()