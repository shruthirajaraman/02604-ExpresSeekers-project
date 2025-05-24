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

# ----------------------- 1. Training & Evaluation -----------------------

# This function takes a model name, model instance, training and testing data, label encoder, and optional feature names.
# It evaluates the model using cross-validation and computes various metrics.
def evaluate_model_all_genes(name, model, X_train, y_train, X_test, y_test, label_encoder, feature_names=None):
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
        # Feature importance plot
    if hasattr(model, 'feature_importances_') and feature_names is not None:
        importances = model.feature_importances_
        indices = np.argsort(importances)[::-1][:3000]  # Get top 3000
        top_features = np.array(feature_names)[indices].tolist()
        
    elif hasattr(model, 'coef_') and feature_names is not None:
        importances = np.abs(model.coef_).sum(axis=0)
        indices = np.argsort(importances)[::-1][:3000]  # Get top 3000
        top_features = np.array(feature_names)[indices].tolist()

    return results, top_features 

# ----------------------- 2. Model Training -----------------------

# This function takes the training and testing data, label encoder, gene set label, and feature names.
# It trains and evaluates multiple models (SVM, Logistic Regression, XGBoost) using the specified gene set.
def train_and_evaluate_models_all_genes(X_train, y_train, X_test, y_test, label_encoder, gene_set_label, feature_names):
    models = {
        'SVM': SVC(kernel='linear', probability=False, cache_size=1000, tol=0.001, max_iter=1000, shrinking=True),
        'LogisticRegression': LogisticRegression(max_iter=1000),
        'XGBoost': XGBClassifier(tree_method='hist', n_jobs=-1,  n_estimators=200, learning_rate=0.1, max_depth=6, subsample=0.8, colsample_bytree=0.8, gamma=0.1, min_child_weight=5, eval_metric='mlogloss')
    }
    all_results = {}
    model_top_genes = {}
    for name, model in models.items():
        print(f"Training {name} with {gene_set_label}...")
        results, top_genes = evaluate_model_all_genes(name, model, X_train, y_train, X_test, y_test, label_encoder, feature_names)
        all_results[f"{name}_{gene_set_label}"] = results['Test']
        model_top_genes[f"{name}_{gene_set_label}"] = top_genes

    return all_results, model_top_genes