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
from rnanorm import CPM 

# ----------------------- 1. Feature Selection -----------------------

# This function takes a Dataframe 
# It sekects akk of the genes
def select_all_genes(df):
    all_genes = set(df.columns.drop('cases.disease_type'))
    available_genes = list(all_genes)
    return available_genes

# ----------------------- 2. Training & Evaluation -----------------------

# This function takes a model name, model instance, training and testing data, label encoder, and optional feature names.
# It evaluates the model using cross-validation and computes various metrics.
def evaluate_model_all_vs_sig(name, model, X_train, y_train, X_test, y_test, label_encoder):
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

    # Print metrics
    print(f"\n{name} - Metrics:")
    print(f"CV Accuracy: {results['CV']['Accuracy']:.4f}")
    print(f"Test Accuracy: {results['Test']['Accuracy']:.4f}")

    return results
    
# ----------------------- 3. Comparison Plot -----------------------

# This function takes two sets of results (for significant genes and all genes) and plots a comparison of the metrics.
def plot_comparison_random_all_vs_sig(results_sig, results_all):
    metrics = ['Accuracy', 'Precision', 'Recall', 'F1', 'MCC']
    df = []
    
    # Process XGBoost (Significant Genes) results
    for metric in metrics:
        df.append({
            'Model': 'XGBoost (Significant Genes)',
            'Metric': metric,
            'Score': results_sig['Test'][metric]
        })
    
    # Process SVM (All Genes) results
    for metric in metrics:
        df.append({
            'Model': 'SVM (All Genes)',
            'Metric': metric,
            'Score': results_all['Test'][metric]
        })
    
    df = pd.DataFrame(df)
    
    plt.figure(figsize=(10, 6))
    sns.barplot(x='Metric', y='Score', hue='Model', data=df)
    plt.title('Model Performance Comparison', fontsize=14)
    plt.ylim(0.8, 1.0)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()