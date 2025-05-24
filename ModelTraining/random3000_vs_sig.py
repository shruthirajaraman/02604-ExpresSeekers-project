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

# This function takes a DataFrame and a list of significant genes, and number of random genes to collect.
# It returns a list of non-significant genes that are present in the DataFrame and do not overlap with the significant genes.
def select_non_overlapping_random_genes(df, sig_genes, num_genes, seed=42):
    all_genes = set(df.columns.drop('cases.disease_type'))
    available_genes = list(all_genes - set(sig_genes))

    if len(available_genes) < num_genes:
        raise ValueError(f"Only {len(available_genes)} non-overlapping genes available")

    random.seed(seed)
    return random.sample(available_genes, num_genes)

# ----------------------- 2. Training & Evaluation -----------------------

# This function takes a model name, model instance, training and testing data, label encoder, and optional feature names.
# It evaluates the model using cross-validation and computes various metrics.
def evaluate_model_random_vs_sig(name, model, X_train, y_train, X_test, y_test, label_encoder):
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

# This function takes the results of two models (random and significant genes) and plots a comparison of their metrics.
def plot_comparison_random_vs_sig(results_random, results_significant):
    metrics = ['Accuracy', 'Precision', 'Recall', 'F1', 'MCC']
    df = []

    for model_label, metrics_dict in {'Random': results_random, 'Significant': results_significant}.items():
        for metric in metrics:
            df.append({
                'Model': model_label,
                'Metric': metric,
                'Score': metrics_dict[metric]
            })

    df = pd.DataFrame(df)
    plt.figure(figsize=(8, 4))
    sns.barplot(data=df, x='Metric', y='Score', hue='Model')
    plt.title('Comparison of Test Metrics: Significant vs Random Genes')
    plt.xticks(rotation=45)
    plt.ylim(0.60, 1.00)
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.show()