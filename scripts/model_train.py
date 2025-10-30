"""
    OVERVIEW

    model_train.py

    Functional Overview:
        This module is responsible for training random forest classifier models
        used for classifying candidate regions in the PICI identification pipeline.
        Contains two independent model training functions optimized for different feature sets.

    Data Processing Pipeline:
        - Load training datasets from database files
        - Data standardization and column name configuration
        - Dataset splitting (70% training set, 30% test set)
        - Model training and feature importance analysis

    Model Parameters:
        - Number of trees: 10,000
        - Random seed: 0 (ensures reproducible results)
        - Parallel computing: Utilizes all available CPU cores
"""

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split

def forest_tran(filepath):
    """
    Load and return the pre-trained random forest classifier model for PICI region feature classification.
    This model utilizes multidimensional features including gene density, hypothetical protein ratio,
    and directional consistency.

    Parameter:
        filepath : str

    Return: sklearn.ensemble.RandomForestClassifier
    """

    # Load the dataset
    url1 = pd.read_table(filepath + '/databases/new_PICI_feature_forest.txt', header=None)
    url1 = pd.DataFrame(url1)
    url1.columns = ['Class_label', 'gene_density', 'gene_hp', 'gene_dis', 'gene_mid']
    Class_label = np.unique(url1['Class_label'])

    info_url = url1.info()
    # Split the training set and test set
    x, y = url1.iloc[:, 1:].values, url1.iloc[:, 0].values
    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.3, random_state=0)
    feat_labels = url1.columns[1:]

    # Train model
    forest = RandomForestClassifier(n_estimators=10000, random_state=0, n_jobs=-1)
    forest.fit(x_train, y_train)
    importances = forest.feature_importances_
    x_columns = url1.columns[1:]
    indices = np.argsort(importances)[::-1]

    # Output important decision features
    for f in range(x_train.shape[1]):
        print("%2d) %-*s %f" % (f + 1, 30, feat_labels[indices[f]], importances[indices[f]]))
    return forest

def forest_tran1(filepath):
    """
    Load the pre-trained random forest model for PICI classification using protein ratios.

    Parameter:
       filepath : str

    Return: sklearn.ensemble.RandomForestClassifier
    """

    url1 = pd.read_table(filepath + '/databases/new_PICI_feature_forest1.txt', header=None)
    url1 = pd.DataFrame(url1)
    url1.columns = ['Class_label', 'gene_hp']
    Class_label = np.unique(url1['Class_label'])

    info_url = url1.info()
    x, y = url1.iloc[:, 1:].values, url1.iloc[:, 0].values
    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.3, random_state=0)
    feat_labels = url1.columns[1:]

    forest1 = RandomForestClassifier(n_estimators=10000, random_state=0, n_jobs=-1)
    forest1.fit(x_train, y_train)
    importances = forest1.feature_importances_
    x_columns = url1.columns[1:]
    indices = np.argsort(importances)[::-1]

    for f in range(x_train.shape[1]):
        print("%2d) %-*s %f" % (f + 1, 30, feat_labels[indices[f]], importances[indices[f]]))
    return forest1
