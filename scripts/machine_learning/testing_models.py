# !/usr/bin/env python
"""
Author: Evy Kuenen
Date: 26-2-2025
Functionality: This script gets the labels for training if a sample was contaminated or not.
                random forest, gradient boost, xgboost, k-nearest neighbours, naive bayes, complement naive bayes
                ,linear discriminant analysis, quadratic discriminant and ensemble classifier.
                This script also plots the random forest and the ensemble, a confusion matrix, the accuracy of all
                models,the f1 score of all models and de mcc score of all models.
Usage:
    Required directory structure:
        train, test and validation file needs to be in test_train_data_ml directory
    Required files:
        train_file = "..\\..\\output_from_scripts\\test_train_data_ml\\train_set.csv"
        test_file = "..\\..\\output_from_scripts\\test_train_data_ml\\test_set.csv"
        val_file = "..\\..\\output_from_scripts\\test_train_data_ml\\val_set.csv"
    Calling script:
        python testing_models.py
"""
import numpy as np
import pandas
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis
from sklearn.ensemble import (RandomForestClassifier, ExtraTreesClassifier, VotingClassifier,
                              GradientBoostingClassifier, BaggingClassifier)
import matplotlib.pyplot as plt
from sklearn.naive_bayes import GaussianNB, ComplementNB
from sklearn.neighbors import KNeighborsClassifier
import xgboost as xgb
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.inspection import permutation_importance
from sklearn.metrics import accuracy_score, f1_score, matthews_corrcoef, confusion_matrix, ConfusionMatrixDisplay, \
    precision_score, classification_report, roc_auc_score
import joblib
from sklearn.tree import plot_tree
from imblearn.over_sampling import SMOTE, SMOTEN
from sklearn.metrics import roc_curve
from sklearn.utils import compute_class_weight


def read_file(train_file, test_file, val_file):
    """
    Reads the files with the training and test data that is going to be used for the models.
    Returns:
        test_data (pandas dataframe): pandas dataframe with the data from the test set.
        train_data (pandas dataframe): pandas dataframe with the data from the training set.
    """
    train_data = pandas.read_csv(train_file, sep=",")
    test_data = pandas.read_csv(test_file, sep=",")
    val_data = pandas.read_csv(val_file, sep=",")
    return test_data, train_data, val_data


def get_labels(train_data, test_data, val_data):
    """
    Gets label if contaminated label is 1 and not contaminated is 0.

    Args:
        train_data (pandas dataframe) : pandas dataframe with the data from the training set.
        test_data (pandas dataframe):  pandas dataframe with the data from the test set.
    Returns:
        train_x_labels (pandas dataframe): pandas dataframe with the rows which will be used to train the model.
        train_y_labels (pandas dataframe): pandas dataframe with the contamination for each row in the train_x_label
        dataframe.
        test_x_labels (pandas dataframe): pandas dataframe with the rows which will be used to test the model.
        test_y_labels (pandas dataframe): pandas dataframe with the contamination for each row in the test_x_label
        dataframe.
        sample,host,virus,ani,correlation,fragmentation,contamination
    """
    column = "contamination"  # label that prediction should have

    # Splits the response(contamination) and the response(the other columns).
    train_x_labels = train_data.loc[:, train_data.columns != column]
    train_y_labels = train_data[column]
    test_x_labels = test_data.loc[:, test_data.columns != column]
    test_y_labels = test_data[column]
    val_x_labels = val_data.loc[:, val_data.columns != column]
    val_y_labels = val_data[column]
    return train_x_labels, train_y_labels, test_x_labels, test_y_labels, val_x_labels, val_y_labels


def train_model_random_forest(train_x_labels, train_y_labels, test_x_labels, test_y_labels, val_x_labels,
                                      val_y_labels):
    """
    Trains and tests the random forest model. Validation is done with k-fold cross validation.
    Also calculates the accuracy, f1 score and MCC score of the model and plots a confusion matrix.
    :param train_x_labels: pandas dataframe with the Feature set for training.
    :param train_y_labels: pandas dataframe with the Labels for training data.
    :param test_x_labels: pandas dataframe with the Feature set for testing.
    :param test_y_labels: pandas dataframe with the Labels for testing data.
    :param val_x_labels: pandas dataframe with the rows which will be used to test the model.
    :param val_y_labels: pandas dataframe with the label of contamination for each row in the test_x_label
         dataframe.
    :return:
        accuracy_score_random_forest (float): shows how many contamination the model predicted correctly.
        randomforest_tree_f1 (float): a combination of precision and recall which is used for evaluating models.
        mcc_rf (float): a score between -1 and 1 which is used for evaluating models.
    """

    print("random forest")
    labels = [0, 1]

    # Feature selectie
    train_x_labels = train_x_labels.iloc[:, -5:].values.astype(np.float32)
    test_x_labels = test_x_labels.iloc[:, -5:].values.astype(np.float32)
    val_x_labels = val_x_labels.iloc[:, -5:].values.astype(np.float32)

    # Labels
    train_y_labels = train_y_labels.values.astype(np.float32)
    test_y_labels = test_y_labels.values.astype(np.float32)
    val_y_labels = val_y_labels.values.astype(np.float32)

    # oversample
    smote = SMOTEN(sampling_strategy=0.25, random_state=42)  # 25% of dataset is minority class
    train_x_resampled, train_y_resampled = smote.fit_resample(train_x_labels, train_y_labels)

    # scale data
    scaler = StandardScaler()
    train_x_scaled = scaler.fit_transform(train_x_resampled)
    test_x_scaled = scaler.transform(test_x_labels)
    val_x_labels = scaler.transform(val_x_labels)

    # give weights to labels
    class_weights = compute_class_weight("balanced", classes=np.array([0, 1]), y=train_y_resampled)
    class_weight_dict = {0: class_weights[0], 1: class_weights[1]}

    # model
    rf = RandomForestClassifier(n_estimators=10000,
                                max_features='sqrt',
                                class_weight=class_weight_dict,
                                bootstrap=True,
                                min_samples_split=5,
                                min_samples_leaf=3,
                                random_state=42)

    # train model
    rf.fit(train_x_scaled, train_y_resampled)

    # === VALIDATIE EVALUATIE ===
    val_probs = rf.predict_proba(val_x_labels)[:, 1]
    fpr, tpr, thresholds = roc_curve(val_y_labels, val_probs)
    optimal_idx = np.argmax(tpr - fpr)
    optimal_threshold = thresholds[optimal_idx]
    print(f"Optimale threshold (gebaseerd op validatie): {optimal_threshold:.2f}")

    val_preds = (val_probs > optimal_threshold).astype(int)

    print("=== Validatieset ===")
    print(classification_report(val_y_labels, val_preds, digits=3))
    cm_val = confusion_matrix(val_y_labels, val_preds)
    ConfusionMatrixDisplay(confusion_matrix=cm_val, display_labels=labels).plot()
    plt.title("Confusion Matrix - Validatie - rf")
    plt.show()

    # === TEST EVALUATIE ===
    test_probs = rf.predict_proba(test_x_scaled)[:, 1]
    test_preds = (test_probs > optimal_threshold).astype(int)

    print("=== Testset ===")
    print(classification_report(test_y_labels, test_preds, digits=3))
    cm_test = confusion_matrix(test_y_labels, test_preds)
    ConfusionMatrixDisplay(confusion_matrix=cm_test, display_labels=labels).plot()
    plt.title("Confusion Matrix - Test - rf")
    plt.show()

    # VALIDATIE: MCC en F1
    val_f1 = f1_score(val_y_labels, val_preds, average='weighted')
    val_mcc = matthews_corrcoef(val_y_labels, val_preds)
    print(f"Validation F1-score (weighted): {val_f1:.3f}")
    print(f"Validation MCC: {val_mcc:.3f}")

    # TEST: MCC en F1
    test_f1 = f1_score(test_y_labels, test_preds, average='weighted')
    test_mcc = matthews_corrcoef(test_y_labels, test_preds)
    print(f"Test F1-score (weighted): {test_f1:.3f}")
    print(f"Test MCC: {test_mcc:.3f}")

    # crossvalidation
    k_folds = StratifiedKFold(n_splits=10, shuffle=True)
    scores_cross_validation = cross_val_score(rf, train_x_scaled, train_y_resampled, cv=k_folds)

    print("Random forest tree")
    print("Average CV Score: ", scores_cross_validation.mean())
    print("Number of CV Scores used in Average: ", len(scores_cross_validation))
    print("")

    # Plots the score of every k-fold.
    plt.plot(range(len(scores_cross_validation)), scores_cross_validation)
    plt.xlabel("K-folds")
    plt.ylabel("Cross validation score")
    plt.title("Cross Validation scores_cross_validation per k-fold")
    plt.show()

    # ROC AUC
    roc_auc_val = roc_auc_score(val_y_labels, val_probs)
    roc_auc_test = roc_auc_score(test_y_labels, test_probs)

    fpr_val, tpr_val, _ = roc_curve(val_y_labels, val_probs)
    fpr_test, tpr_test, _ = roc_curve(test_y_labels, test_probs)

    plt.figure(figsize=(8, 6))
    plt.plot(fpr_val, tpr_val, label=f'Validatie ROC curve (AUC = {roc_auc_val:.2f})')
    plt.plot(fpr_test, tpr_test, label=f'Test ROC curve (AUC = {roc_auc_test:.2f})', linestyle='--')
    plt.plot([0, 1], [0, 1], color='gray', linestyle=':', label='Random Classifier')

    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve - rf')
    plt.legend(loc='lower right')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    print(f"ROC AUC (validatie): {roc_auc_val:.3f}")
    print(f"ROC AUC (test): {roc_auc_test:.3f}")

    # Calculates the importance of every feature used to come to the predictions.
    importance = rf.feature_importances_
    features = ['ani', "correlation", 'fragmentation', 'average contig length', 'contig_count']
    # Summarize feature importance.
    for i, v in enumerate(importance):
        print('Feature:  %0d, Score: %.5f' % (i, v))
    # Plot feature importance.
    plt.bar([x for x in features], importance)
    plt.xlabel("Features")
    plt.ylabel("Importance")
    plt.title("Feature Importance for the random forest model")
    plt.show()

    plt.figure(figsize=(15, 10))
    plot_tree(rf.estimators_[0], feature_names=features, filled=True, rounded=True, fontsize=8)
    plt.show()

    randomforest_tree_f1 = f1_score(test_y_labels, test_preds, average="weighted")
    accuracy_score_random_forest = accuracy_score(test_y_labels, test_preds)
    mcc_rf = matthews_corrcoef(test_y_labels, test_preds)

    # # Saves the model into a file.
    # joblib.dump(rf, "rf_2xaugment_mcc0.37.joblib")

    return accuracy_score_random_forest, randomforest_tree_f1, mcc_rf


def train_model_mlp(train_x_labels, train_y_labels, test_x_labels, test_y_labels, val_x_labels, val_y_labels):
    """
    Trains and tests a Multi-Layer Perceptron (MLP) model with class imbalance handling.
    Calculates accuracy, F1 score, MCC score, and plots a confusion matrix.
    :param train_x_labels: pandas dataframe with the Feature set for training.
    :param train_y_labels: pandas dataframe with the Labels for training data.
    :param test_x_labels: pandas dataframe with the Feature set for testing.
    :param test_y_labels: pandas dataframe with the Labels for testing data.
    :param val_x_labels: pandas dataframe with the rows which will be used to test the model.
    :param val_y_labels: pandas dataframe with the label of contamination for each row in the test_x_label
         dataframe.
    :return:
        accuracy_score_mlp (float): Accuracy of the model.
        mlp_f1 (float): F1-score (weighted) for model evaluation.
        mcc_mlp (float): Matthews correlation coefficient.
    """

    labels = [0, 1]
    print("Multi-Layer Perceptron (MLP)")

    # Feature selectie
    train_x_labels = train_x_labels.iloc[:, -5:].values.astype(np.float32)
    test_x_labels = test_x_labels.iloc[:, -5:].values.astype(np.float32)
    val_x_labels = val_x_labels.iloc[:, -5:].values.astype(np.float32)

    # Labels
    train_y_labels = train_y_labels.values.astype(np.float32)
    test_y_labels = test_y_labels.values.astype(np.float32)
    val_y_labels = val_y_labels.values.astype(np.float32)

    # oversample
    smote = SMOTEN(sampling_strategy=0.1, random_state=42)
    train_x_resampled, train_y_resampled = smote.fit_resample(train_x_labels, train_y_labels)
    train_y_int = train_y_resampled.astype(int)

    # normalize
    scaler = StandardScaler()
    train_x_scaled = scaler.fit_transform(train_x_resampled)
    test_x_scaled = scaler.transform(test_x_labels)
    val_x_labels = scaler.transform(val_x_labels)

    # give weights to labels
    class_weights = compute_class_weight(class_weight='balanced', classes=np.array([0, 1]), y=train_y_int)
    pos_weight_value = class_weights[1] / class_weights[0]
    print(f"Berekenede pos_weight: {pos_weight_value:.2f}")

    # train model
    model = MLPClassifier(random_state=42,
                          validation_fraction=0.3,
                          learning_rate_init=0.001,
                          )

    # train model
    model.fit(train_x_scaled, train_y_resampled)

    # === VALIDATIE EVALUATIE ===
    val_probs = model.predict_proba(val_x_labels)[:, 1]
    fpr, tpr, thresholds = roc_curve(val_y_labels, val_probs)
    optimal_idx = np.argmax(tpr - fpr)
    optimal_threshold = thresholds[optimal_idx]
    print(f"Optimale threshold (gebaseerd op validatie): {optimal_threshold:.2f}")

    val_preds = (val_probs > optimal_threshold).astype(int)

    print("=== Validatieset ===")
    print(classification_report(val_y_labels, val_preds, digits=3))
    cm_val = confusion_matrix(val_y_labels, val_preds)
    ConfusionMatrixDisplay(confusion_matrix=cm_val, display_labels=labels).plot()
    plt.title("Confusion Matrix - Validatie - mlp")
    plt.show()

    # === TEST EVALUATIE ===
    test_probs = model.predict_proba(test_x_scaled)[:, 1]
    test_preds = (test_probs > optimal_threshold).astype(int)

    print("=== Testset ===")
    print(classification_report(test_y_labels, test_preds, digits=3))
    cm_test = confusion_matrix(test_y_labels, test_preds)
    ConfusionMatrixDisplay(confusion_matrix=cm_test, display_labels=labels).plot()
    plt.title("Confusion Matrix - Test - mlp")
    plt.show()

    # VALIDATIE: MCC en F1
    val_f1 = f1_score(val_y_labels, val_preds, average='weighted')
    val_mcc = matthews_corrcoef(val_y_labels, val_preds)
    print(f"Validation F1-score (weighted): {val_f1:.3f}")
    print(f"Validation MCC: {val_mcc:.3f}")

    # TEST: MCC en F1
    test_f1 = f1_score(test_y_labels, test_preds, average='weighted')
    test_mcc = matthews_corrcoef(test_y_labels, test_preds)
    print(f"Test F1-score (weighted): {test_f1:.3f}")
    print(f"Test MCC: {test_mcc:.3f}")

    # ROC AUC
    roc_auc_val = roc_auc_score(val_y_labels, val_probs)
    roc_auc_test = roc_auc_score(test_y_labels, test_probs)

    fpr_val, tpr_val, _ = roc_curve(val_y_labels, val_probs)
    fpr_test, tpr_test, _ = roc_curve(test_y_labels, test_probs)

    plt.figure(figsize=(8, 6))
    plt.plot(fpr_val, tpr_val, label=f'Validatie ROC curve (AUC = {roc_auc_val:.2f})')
    plt.plot(fpr_test, tpr_test, label=f'Test ROC curve (AUC = {roc_auc_test:.2f})', linestyle='--')
    plt.plot([0, 1], [0, 1], color='gray', linestyle=':', label='Random Classifier')

    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve - mlp')
    plt.legend(loc='lower right')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    print(f"ROC AUC (validatie): {roc_auc_val:.3f}")
    print(f"ROC AUC (test): {roc_auc_test:.3f}")

    # crossvalidation
    k_folds = StratifiedKFold(n_splits=10, shuffle=True)
    scores_cross_validation = cross_val_score(model, train_x_scaled, train_y_resampled, cv=k_folds)

    print("Random forest tree")
    print("Average CV Score: ", scores_cross_validation.mean())
    print("Number of CV Scores used in Average: ", len(scores_cross_validation))
    print("")

    mlp_f1 = f1_score(test_y_labels, test_preds, average="weighted")
    accuracy_score_mlp = accuracy_score(test_y_labels, test_preds)
    mcc_mlp = matthews_corrcoef(test_y_labels, test_preds)

    return accuracy_score_mlp, mlp_f1, mcc_mlp


def train_model_naive_bayes(train_x_labels, train_y_labels, test_x_labels, test_y_labels, val_x_labels, val_y_labels):
    """
    Trains and tests the Naive Bayes model. Also calculates the accuracy, f1 score
    and MCC score of the model and plots a confusion matrix.
    :param train_x_labels: pandas dataframe with the Feature set for training.
    :param train_y_labels: pandas dataframe with the Labels for training data.
    :param test_x_labels: pandas dataframe with the Feature set for testing.
    :param test_y_labels: pandas dataframe with the Labels for testing data.
    :param val_x_labels: pandas dataframe with the rows which will be used to test the model.
    :param val_y_labels: pandas dataframe with the label of contamination for each row in the test_x_label
         dataframe.
    :return:
        accuracy_score_naive_bayes (float): shows how many contamination the model predicted correctly.
        naive_bayes_f1 (float): a combination of precision and recall which is used for evaluating models.
        mcc_naive (float): a score between -1 and 1 which is used for evaluating models.
    """
    print("Naive bayes")
    labels = [0, 1]

    # Feature selectie
    train_x_labels = train_x_labels.iloc[:, -5:].values.astype(np.float32)
    test_x_labels = test_x_labels.iloc[:, -5:].values.astype(np.float32)
    val_x_labels = val_x_labels.iloc[:, -5:].values.astype(np.float32)

    # Labels
    train_y_labels = train_y_labels.astype(int)
    test_y_labels = test_y_labels.astype(int)
    val_y_labels = val_y_labels.astype(int)

    # oversample minority class
    smote = SMOTEN(sampling_strategy=0.15, random_state=42)  # 15% of dataset is minority class
    train_x_resampled, train_y_resampled = smote.fit_resample(train_x_labels, train_y_labels)

    # scale
    scaler = MinMaxScaler()
    train_x_resampled = scaler.fit_transform(train_x_resampled)
    test_x_labels = scaler.transform(test_x_labels)
    val_x_labels = scaler.transform(val_x_labels)

    # train model
    classifier = GaussianNB(priors=None, var_smoothing=1e-9)
    classifier.fit(train_x_resampled, train_y_resampled)

    # === VALIDATIE EVALUATIE ===
    val_probs = classifier.predict_proba(val_x_labels)[:, 1]
    fpr, tpr, thresholds = roc_curve(val_y_labels, val_probs)
    optimal_idx = np.argmax(tpr - fpr)
    optimal_threshold = thresholds[optimal_idx]
    print(f"Optimale threshold (gebaseerd op validatie): {optimal_threshold:.2f}")

    val_preds = (val_probs > optimal_threshold).astype(int)

    print("=== Validatieset ===")
    print(classification_report(val_y_labels, val_preds, digits=3))
    cm_val = confusion_matrix(val_y_labels, val_preds)
    ConfusionMatrixDisplay(confusion_matrix=cm_val, display_labels=labels).plot()
    plt.title("Confusion Matrix - Validatie - nb")
    plt.show()

    # === TEST EVALUATIE ===
    test_probs = classifier.predict_proba(test_x_labels)[:, 1]
    test_preds = (test_probs > optimal_threshold).astype(int)

    print("=== Testset ===")
    print(classification_report(test_y_labels, test_preds, digits=3))
    cm_test = confusion_matrix(test_y_labels, test_preds)
    ConfusionMatrixDisplay(confusion_matrix=cm_test, display_labels=labels).plot()
    plt.title("Confusion Matrix - Test - nb")
    plt.show()

    # VALIDATIE: MCC en F1
    val_f1 = f1_score(val_y_labels, val_preds, average='weighted')
    val_mcc = matthews_corrcoef(val_y_labels, val_preds)
    print(f"Validation F1-score (weighted): {val_f1:.3f}")
    print(f"Validation MCC: {val_mcc:.3f}")

    # TEST: MCC en F1
    test_f1 = f1_score(test_y_labels, test_preds, average='weighted')
    test_mcc = matthews_corrcoef(test_y_labels, test_preds)
    print(f"Test F1-score (weighted): {test_f1:.3f}")
    print(f"Test MCC: {test_mcc:.3f}")

    naive_bayes_f1 = f1_score(test_y_labels, test_preds, average="weighted")
    accuracy_score_naive_bayes = accuracy_score(test_y_labels, test_preds)
    mcc_naive = matthews_corrcoef(test_y_labels, test_preds)
    print("Accuracy score: ", (accuracy_score_naive_bayes / 1) * 100, "%")
    print("F1 score: ", naive_bayes_f1)
    print("MCC score: ", mcc_naive)
    print("precision: ", precision_score(test_y_labels, test_preds))
    print("")

    # ROC AUC
    roc_auc_val = roc_auc_score(val_y_labels, val_probs)
    roc_auc_test = roc_auc_score(test_y_labels, test_probs)

    fpr_val, tpr_val, _ = roc_curve(val_y_labels, val_probs)
    fpr_test, tpr_test, _ = roc_curve(test_y_labels, test_probs)

    plt.figure(figsize=(8, 6))
    plt.plot(fpr_val, tpr_val, label=f'Validatie ROC curve (AUC = {roc_auc_val:.2f})')
    plt.plot(fpr_test, tpr_test, label=f'Test ROC curve (AUC = {roc_auc_test:.2f})', linestyle='--')
    plt.plot([0, 1], [0, 1], color='gray', linestyle=':', label='Random Classifier')

    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve - nb')
    plt.legend(loc='lower right')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    print(f"ROC AUC (validatie): {roc_auc_val:.3f}")
    print(f"ROC AUC (test): {roc_auc_test:.3f}")

    # crossvalidation
    k_folds = StratifiedKFold(n_splits=10, shuffle=True)
    scores_cross_validation = cross_val_score(classifier, train_x_resampled, train_y_resampled, cv=k_folds)

    print("Random forest tree")
    print("Average CV Score: ", scores_cross_validation.mean())
    print("Number of CV Scores used in Average: ", len(scores_cross_validation))
    print("")
    return accuracy_score_naive_bayes, naive_bayes_f1, mcc_naive


def train_model_complementnaive_bayes(train_x_labels, train_y_labels, test_x_labels, test_y_labels, val_x_labels,
                                      val_y_labels):
    """
    Trains and tests the complement Naive Bayes model. Also calculates the accuracy, f1 score
     and MCC score of the model and plots a confusion matrix.
    :param train_x_labels: pandas dataframe with the Feature set for training.
    :param train_y_labels: pandas dataframe with the Labels for training data.
    :param test_x_labels: pandas dataframe with the Feature set for testing.
    :param test_y_labels: pandas dataframe with the Labels for testing data.
    :param val_x_labels: pandas dataframe with the rows which will be used to test the model.
    :param val_y_labels: pandas dataframe with the label of contamination for each row in the test_x_label
         dataframe.
    :return:
        accuracy_score_complementnaive_bayes (float): shows how many contamination the model predicted correctly.
         complementnaive_bayes_f1 (float): a combination of precision and recall which is used for evaluating models.
         mcc_complementnaive (float): a score between -1 and 1 which is used for evaluating models.
    """
    print('Complement Naive Bayes model')
    labels = [0, 1]

    # Feature selectie
    train_x_labels = train_x_labels.iloc[:, -5:].values.astype(np.float32)
    test_x_labels = test_x_labels.iloc[:, -5:].values.astype(np.float32)
    val_x_labels = val_x_labels.iloc[:, -5:].values.astype(np.float32)

    # Labels
    train_y_labels = train_y_labels.astype(int)
    test_y_labels = test_y_labels.astype(int)
    val_y_labels = val_y_labels.astype(int)

    # Oversample minority class in train set
    smote = SMOTE(sampling_strategy=0.15, random_state=42)
    train_x_resampled, train_y_resampled = smote.fit_resample(train_x_labels, train_y_labels)

    # Scaling
    scaler = MinMaxScaler()
    train_x_resampled = scaler.fit_transform(train_x_resampled)
    test_x_labels = scaler.transform(test_x_labels)
    val_x_labels = scaler.transform(val_x_labels)

    # Train model
    classifier = ComplementNB(fit_prior=True, norm=True)
    classifier.fit(train_x_resampled, train_y_resampled)

    # === VALIDATIE EVALUATIE ===
    val_probs = classifier.predict_proba(val_x_labels)[:, 1]
    fpr, tpr, thresholds = roc_curve(val_y_labels, val_probs)
    optimal_idx = np.argmax(tpr - fpr)
    optimal_threshold = thresholds[optimal_idx]
    print(f"Optimale threshold (gebaseerd op validatie): {optimal_threshold:.2f}")

    val_preds = (val_probs > optimal_threshold).astype(int)

    print("=== Validatieset ===")
    print(classification_report(val_y_labels, val_preds, digits=3))
    cm_val = confusion_matrix(val_y_labels, val_preds)
    ConfusionMatrixDisplay(confusion_matrix=cm_val, display_labels=labels).plot()
    plt.title("Confusion Matrix - Validatie - cnb")
    plt.show()

    # === TEST EVALUATIE ===
    test_probs = classifier.predict_proba(test_x_labels)[:, 1]
    test_preds = (test_probs > optimal_threshold).astype(int)

    print("=== Testset ===")
    print(classification_report(test_y_labels, test_preds, digits=3))
    cm_test = confusion_matrix(test_y_labels, test_preds)
    ConfusionMatrixDisplay(confusion_matrix=cm_test, display_labels=labels).plot()
    plt.title("Confusion Matrix - Test - cnb")
    plt.show()

    # VALIDATIE: MCC en F1
    val_f1 = f1_score(val_y_labels, val_preds, average='weighted')
    val_mcc = matthews_corrcoef(val_y_labels, val_preds)
    print(f"Validation F1-score (weighted): {val_f1:.3f}")
    print(f"Validation MCC: {val_mcc:.3f}")

    # TEST: MCC en F1
    test_f1 = f1_score(test_y_labels, test_preds, average='weighted')
    test_mcc = matthews_corrcoef(test_y_labels, test_preds)
    print(f"Test F1-score (weighted): {test_f1:.3f}")
    print(f"Test MCC: {test_mcc:.3f}")

    # Return kernscores
    test_accuracy = accuracy_score(test_y_labels, test_preds)
    test_f1 = f1_score(test_y_labels, test_preds, average="weighted")
    test_mcc = matthews_corrcoef(test_y_labels, test_preds)

    # ROC AUC
    roc_auc_val = roc_auc_score(val_y_labels, val_probs)
    roc_auc_test = roc_auc_score(test_y_labels, test_probs)

    fpr_val, tpr_val, _ = roc_curve(val_y_labels, val_probs)
    fpr_test, tpr_test, _ = roc_curve(test_y_labels, test_probs)

    plt.figure(figsize=(8, 6))
    plt.plot(fpr_val, tpr_val, label=f'Validatie ROC curve (AUC = {roc_auc_val:.2f})')
    plt.plot(fpr_test, tpr_test, label=f'Test ROC curve (AUC = {roc_auc_test:.2f})', linestyle='--')
    plt.plot([0, 1], [0, 1], color='gray', linestyle=':', label='Random Classifier')

    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve cnb')
    plt.legend(loc='lower right')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    print(f"ROC AUC (validatie): {roc_auc_val:.3f}")
    print(f"ROC AUC (test): {roc_auc_test:.3f}")

    # crossvalidation
    k_folds = StratifiedKFold(n_splits=10, shuffle=True)
    scores_cross_validation = cross_val_score(classifier, train_x_resampled, train_y_resampled, cv=k_folds)

    print("Random forest tree")
    print("Average CV Score: ", scores_cross_validation.mean())
    print("Number of CV Scores used in Average: ", len(scores_cross_validation))
    print("")

    # feature importance
    result = permutation_importance(classifier, val_x_labels, val_y_labels, n_repeats=10, random_state=42)
    sorted_idx = result.importances_mean.argsort()[::-1]
    feature_names = ['ani', "correlation", 'fragmentation', 'average contig length', "contig_count"]
    for i in sorted_idx[:10]:
        print(f"{feature_names[i]}: {result.importances_mean[i]:.4f}")
    return test_accuracy, test_f1, test_mcc


def linear_discriminant_analysis_model(train_x_labels, train_y_labels, test_x_labels, test_y_labels, val_x_labels,
                                      val_y_labels):
    """
    Trains and tests the Linear Discriminant model. Also calculates the accuracy, f1 score
    and MCC score of the model and plots a confusion matrix
    :param train_x_labels: pandas dataframe with the Feature set for training.
    :param train_y_labels: pandas dataframe with the Labels for training data.
    :param test_x_labels: pandas dataframe with the Feature set for testing.
    :param test_y_labels: pandas dataframe with the Labels for testing data.
    :param val_x_labels: pandas dataframe with the rows which will be used to test the model.
    :param val_y_labels: pandas dataframe with the label of contamination for each row in the test_x_label
            dataframe.
    :return:
        accuracy_score_lda (float): shows how many contamination the model predicted correctly.
        f1_lda (float): a combination of precision and recall which is used for evaluating models.
        mcc_lda (float): a score between -1 and 1 which is used for evaluating models.
    """
    print("linear discriminant")
    labels = [0, 1]

    # Feature selectie
    train_x_labels = train_x_labels.iloc[:, -5:].values.astype(np.float32)
    test_x_labels = test_x_labels.iloc[:, -5:].values.astype(np.float32)
    val_x_labels = val_x_labels.iloc[:, -5:].values.astype(np.float32)

    # Labels
    train_y_labels = train_y_labels.astype(int)
    test_y_labels = test_y_labels.astype(int)
    val_y_labels = val_y_labels.astype(int)

    # oversample minority
    smote = SMOTEN(sampling_strategy=0.1, random_state=42)  # 10% of dataset is minority class
    train_x_resampled, train_y_resampled = smote.fit_resample(train_x_labels, train_y_labels)

    # scale
    scaler = MinMaxScaler()
    train_x_resampled = scaler.fit_transform(train_x_resampled)
    test_x_labels = scaler.transform(test_x_labels)
    val_x_labels = scaler.transform(val_x_labels)

    # train model
    lda = LinearDiscriminantAnalysis(solver="lsqr")
    lda.fit(train_x_resampled, train_y_resampled)

    # === VALIDATIE EVALUATIE ===
    val_probs = lda.predict_proba(val_x_labels)[:, 1]
    fpr, tpr, thresholds = roc_curve(val_y_labels, val_probs)
    optimal_idx = np.argmax(tpr - fpr)
    optimal_threshold = thresholds[optimal_idx]
    print(f"Optimale threshold (gebaseerd op validatie): {optimal_threshold:.2f}")

    val_preds = (val_probs > optimal_threshold).astype(int)

    print("=== Validatieset ===")
    print(classification_report(val_y_labels, val_preds, digits=3))
    cm_val = confusion_matrix(val_y_labels, val_preds)
    ConfusionMatrixDisplay(confusion_matrix=cm_val, display_labels=labels).plot()
    plt.title("Confusion Matrix - Validatie - ld")
    plt.show()

    # === TEST EVALUATIE ===
    test_probs = lda.predict_proba(test_x_labels)[:, 1]
    test_preds = (test_probs > optimal_threshold).astype(int)

    print("=== Testset ===")
    print(classification_report(test_y_labels, test_preds, digits=3))
    cm_test = confusion_matrix(test_y_labels, test_preds)
    ConfusionMatrixDisplay(confusion_matrix=cm_test, display_labels=labels).plot()
    plt.title("Confusion Matrix - Test - ld")
    plt.show()

    # VALIDATIE: MCC en F1
    val_f1 = f1_score(val_y_labels, val_preds, average='weighted')
    val_mcc = matthews_corrcoef(val_y_labels, val_preds)
    print(f"Validation F1-score (weighted): {val_f1:.3f}")
    print(f"Validation MCC: {val_mcc:.3f}")

    # TEST: MCC en F1
    test_f1 = f1_score(test_y_labels, test_preds, average='weighted')
    test_mcc = matthews_corrcoef(test_y_labels, test_preds)
    print(f"Test F1-score (weighted): {test_f1:.3f}")
    print(f"Test MCC: {test_mcc:.3f}")

    accuracy_score_lda = accuracy_score(test_y_labels, test_preds)
    f1_lda = f1_score(test_y_labels, test_preds, average="weighted")
    mcc_lda = matthews_corrcoef(test_y_labels, test_preds)
    print("Accuracy score: ", (accuracy_score_lda / 1) * 100, "%")
    print("F1 score: ", f1_lda)
    print("MCC score: ", mcc_lda)
    print("precision: ", precision_score(test_y_labels, test_preds))
    print("")

    # ROC AUC
    roc_auc_val = roc_auc_score(val_y_labels, val_probs)
    roc_auc_test = roc_auc_score(test_y_labels, test_preds)

    fpr_val, tpr_val, _ = roc_curve(val_y_labels, val_probs)
    fpr_test, tpr_test, _ = roc_curve(test_y_labels, test_probs)

    plt.figure(figsize=(8, 6))
    plt.plot(fpr_val, tpr_val, label=f'Validatie ROC curve (AUC = {roc_auc_val:.2f})')
    plt.plot(fpr_test, tpr_test, label=f'Test ROC curve (AUC = {roc_auc_test:.2f})', linestyle='--')
    plt.plot([0, 1], [0, 1], color='gray', linestyle=':', label='Random Classifier')

    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve - lda')
    plt.legend(loc='lower right')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    print(f"ROC AUC (validatie): {roc_auc_val:.3f}")
    print(f"ROC AUC (test): {roc_auc_test:.3f}")

    # crossvalidation
    k_folds = StratifiedKFold(n_splits=10, shuffle=True)
    scores_cross_validation = cross_val_score(lda, train_x_resampled, train_y_resampled, cv=k_folds)

    print("Random forest tree")
    print("Average CV Score: ", scores_cross_validation.mean())
    print("Number of CV Scores used in Average: ", len(scores_cross_validation))
    print("")

    return accuracy_score_lda, f1_lda, mcc_lda


def quadratic_discriminant_model(train_x_labels, train_y_labels, test_x_labels, test_y_labels, val_x_labels,
                                      val_y_labels):
    """
    Trains and tests the Quadratic Discriminant model. Also calculates the accuracy, f1 score
    and MCC score of the model and plots a confusion matrix.
    :param train_x_labels: pandas dataframe with the Feature set for training.
    :param train_y_labels: pandas dataframe with the Labels for training data.
    :param test_x_labels: pandas dataframe with the Feature set for testing.
    :param test_y_labels: pandas dataframe with the Labels for testing data.
    :param val_x_labels: pandas dataframe with the rows which will be used to test the model.
    :param val_y_labels: pandas dataframe with the label of contamination for each row in the test_x_label
            dataframe.
    :return:
        accuracy_score_qda (float): shows how many contamination the model predicted correctly.
        f1_qda (float): a combination of precision and recall which is used for evaluating models.
        mcc_qda (float): a score between -1 and 1 which is used for evaluating models.
    """
    print("quadratic discriminant")
    labels = [0, 1]

    # Feature selectie
    train_x_labels = train_x_labels.iloc[:, -5:].values.astype(np.float32)
    test_x_labels = test_x_labels.iloc[:, -5:].values.astype(np.float32)
    val_x_labels = val_x_labels.iloc[:, -5:].values.astype(np.float32)

    # Labels
    train_y_labels = train_y_labels.astype(int)
    test_y_labels = test_y_labels.astype(int)
    val_y_labels = val_y_labels.astype(int)

    # oversample minority class
    smote = SMOTEN(sampling_strategy=0.1, random_state=42)  # 10% of dataset is minority class
    train_x_resampled, train_y_resampled = smote.fit_resample(train_x_labels, train_y_labels)

    # scale
    scaler = MinMaxScaler()
    train_x_resampled = scaler.fit_transform(train_x_resampled)
    test_x_labels = scaler.transform(test_x_labels)
    val_x_labels = scaler.transform(val_x_labels)

    # train model
    qda = QuadraticDiscriminantAnalysis()
    qda.fit(train_x_resampled, train_y_resampled)

    # === VALIDATIE EVALUATIE ===
    val_probs = qda.predict_proba(val_x_labels)[:, 1]
    fpr, tpr, thresholds = roc_curve(val_y_labels, val_probs)
    optimal_idx = np.argmax(tpr - fpr)
    optimal_threshold = thresholds[optimal_idx]
    print(f"Optimale threshold (gebaseerd op validatie): {optimal_threshold:.2f}")

    val_preds = (val_probs > optimal_threshold).astype(int)

    print("=== Validatieset ===")
    print(classification_report(val_y_labels, val_preds, digits=3))
    cm_val = confusion_matrix(val_y_labels, val_preds)
    ConfusionMatrixDisplay(confusion_matrix=cm_val, display_labels=labels).plot()
    plt.title("Confusion Matrix - Validatie - qd")
    plt.show()

    # === TEST EVALUATIE ===
    test_probs = qda.predict_proba(test_x_labels)[:, 1]
    test_preds = (test_probs > optimal_threshold).astype(int)

    print("=== Testset ===")
    print(classification_report(test_y_labels, test_preds, digits=3))
    cm_test = confusion_matrix(test_y_labels, test_preds)
    ConfusionMatrixDisplay(confusion_matrix=cm_test, display_labels=labels).plot()
    plt.title("Confusion Matrix - Test - qd")
    plt.show()

    # VALIDATIE: MCC en F1
    val_f1 = f1_score(val_y_labels, val_preds, average='weighted')
    val_mcc = matthews_corrcoef(val_y_labels, val_preds)
    print(f"Validation F1-score (weighted): {val_f1:.3f}")
    print(f"Validation MCC: {val_mcc:.3f}")

    # TEST: MCC en F1
    test_f1 = f1_score(test_y_labels, test_preds, average='weighted')
    test_mcc = matthews_corrcoef(test_y_labels, test_preds)
    print(f"Test F1-score (weighted): {test_f1:.3f}")
    print(f"Test MCC: {test_mcc:.3f}")

    accuracy_score_qda = accuracy_score(test_y_labels, test_preds)
    f1_qda = f1_score(test_y_labels, test_preds, average="weighted")
    mcc_qda = matthews_corrcoef(test_y_labels, test_preds)
    print("Accuracy score: ", (accuracy_score_qda / 1) * 100, "%")
    print("F1 score: ", f1_qda)
    print("MCC score: ", mcc_qda)
    print("precision: ", precision_score(test_y_labels, test_preds))
    print("")

    # ROC AUC
    roc_auc_val = roc_auc_score(val_y_labels, val_probs)
    roc_auc_test = roc_auc_score(test_y_labels, test_probs)

    fpr_val, tpr_val, _ = roc_curve(val_y_labels, val_probs)
    fpr_test, tpr_test, _ = roc_curve(test_y_labels, test_probs)

    plt.figure(figsize=(8, 6))
    plt.plot(fpr_val, tpr_val, label=f'Validatie ROC curve (AUC = {roc_auc_val:.2f})')
    plt.plot(fpr_test, tpr_test, label=f'Test ROC curve (AUC = {roc_auc_test:.2f})', linestyle='--')
    plt.plot([0, 1], [0, 1], color='gray', linestyle=':', label='Random Classifier')

    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve - qda')
    plt.legend(loc='lower right')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    print(f"ROC AUC (validatie): {roc_auc_val:.3f}")
    print(f"ROC AUC (test): {roc_auc_test:.3f}")

    # crossvalidation
    k_folds = StratifiedKFold(n_splits=10, shuffle=True)
    scores_cross_validation = cross_val_score(qda, train_x_resampled, train_y_resampled, cv=k_folds)

    print("Random forest tree")
    print("Average CV Score: ", scores_cross_validation.mean())
    print("Number of CV Scores used in Average: ", len(scores_cross_validation))
    print("")
    return accuracy_score_qda, f1_qda, mcc_qda


def gradient_boost(train_x_labels, train_y_labels, test_x_labels, test_y_labels, val_x_labels,
                                      val_y_labels):
    """
    Trains and tests the GradientBoost model. Also calculates the accuracy, f1 score
    and MCC score of the model and plots a confusion matrix.
    :param train_x_labels: pandas dataframe with the Feature set for training.
    :param train_y_labels: pandas dataframe with the Labels for training data.
    :param test_x_labels: pandas dataframe with the Feature set for testing.
    :param test_y_labels: pandas dataframe with the Labels for testing data.
    :param val_x_labels: pandas dataframe with the rows which will be used to test the model.
    :param val_y_labels: pandas dataframe with the label of contamination for each row in the test_x_label
            dataframe.
    :return:
        gb_accuracy (float): shows how much contamination the model predicted correctly.
        gb_f1 (float): a combination of precision and recall which is used for evaluating models.
        mcc_gb (float): a score between -1 and 1 which is used for evaluating models.
    """
    print("Gradient Boost")
    labels = [0, 1]

    # Feature selectie
    train_x_labels = train_x_labels.iloc[:, -5:].values.astype(np.float32)
    test_x_labels = test_x_labels.iloc[:, -5:].values.astype(np.float32)
    val_x_labels = val_x_labels.iloc[:, -5:].values.astype(np.float32)

    # Labels
    train_y_labels = train_y_labels.astype(int)
    test_y_labels = test_y_labels.astype(int)
    val_y_labels = val_y_labels.astype(int)

    # oversample monority class
    smote = SMOTE(sampling_strategy=0.25, random_state=42)  # 25% of dataset is minority class
    train_x_resampled, train_y_resampled = smote.fit_resample(train_x_labels, train_y_labels)

    # scale
    scaler = MinMaxScaler()
    train_x_resampled = scaler.fit_transform(train_x_resampled)
    test_x_labels = scaler.transform(test_x_labels)
    val_x_labels = scaler.transform(val_x_labels)

    # train model
    gb = GradientBoostingClassifier(n_estimators=10000,
                                    max_features=None,
                                    subsample=0.8,
                                    learning_rate=0.01,
                                    )
    gb.fit(train_x_resampled, train_y_resampled)

    # === VALIDATIE EVALUATIE ===
    val_probs = gb.predict_proba(val_x_labels)[:, 1]
    fpr, tpr, thresholds = roc_curve(val_y_labels, val_probs)
    optimal_idx = np.argmax(tpr - fpr)
    optimal_threshold = thresholds[optimal_idx]
    print(f"Optimale threshold (gebaseerd op validatie): {optimal_threshold:.2f}")

    val_preds = (val_probs > optimal_threshold).astype(int)

    print("=== Validatieset ===")
    print(classification_report(val_y_labels, val_preds, digits=3))
    cm_val = confusion_matrix(val_y_labels, val_preds)
    ConfusionMatrixDisplay(confusion_matrix=cm_val, display_labels=labels).plot()
    plt.title("Confusion Matrix - Validatie - gb")
    plt.show()

    # === TEST EVALUATIE ===
    test_probs = gb.predict_proba(test_x_labels)[:, 1]
    test_preds = (test_probs > optimal_threshold).astype(int)

    print("=== Testset ===")
    print(classification_report(test_y_labels, test_preds, digits=3))
    cm_test = confusion_matrix(test_y_labels, test_preds)
    ConfusionMatrixDisplay(confusion_matrix=cm_test, display_labels=labels).plot()
    plt.title("Confusion Matrix - Test - gb")
    plt.show()

    # VALIDATIE: MCC en F1
    val_f1 = f1_score(val_y_labels, val_preds, average='weighted')
    val_mcc = matthews_corrcoef(val_y_labels, val_preds)
    print(f"Validation F1-score (weighted): {val_f1:.3f}")
    print(f"Validation MCC: {val_mcc:.3f}")

    # TEST: MCC en F1
    test_f1 = f1_score(test_y_labels, test_preds, average='weighted')
    test_mcc = matthews_corrcoef(test_y_labels, test_preds)
    print(f"Test F1-score (weighted): {test_f1:.3f}")
    print(f"Test MCC: {test_mcc:.3f}")

    # Calculates the accuracy, f1 and mcc scores and prints them.
    gb_f1 = f1_score(test_y_labels, test_preds, average="weighted")
    gb_accuracy = accuracy_score(test_y_labels, test_preds)
    mcc_gb = matthews_corrcoef(test_y_labels, test_preds)
    print("Accuracy score: ", (gb_accuracy / 1) * 100, "%")
    print("F1 score: ", gb_f1)
    print("mcc: ", mcc_gb)
    print("precision: ", precision_score(test_y_labels, test_preds))
    print("")
    # ROC AUC
    roc_auc_val = roc_auc_score(val_y_labels, val_probs)
    roc_auc_test = roc_auc_score(test_y_labels, test_probs)

    fpr_val, tpr_val, _ = roc_curve(val_y_labels, val_probs)
    fpr_test, tpr_test, _ = roc_curve(test_y_labels, test_probs)

    plt.figure(figsize=(8, 6))
    plt.plot(fpr_val, tpr_val, label=f'Validatie ROC curve (AUC = {roc_auc_val:.2f})')
    plt.plot(fpr_test, tpr_test, label=f'Test ROC curve (AUC = {roc_auc_test:.2f})', linestyle='--')
    plt.plot([0, 1], [0, 1], color='gray', linestyle=':', label='Random Classifier')

    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve - gb')
    plt.legend(loc='lower right')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    print(f"ROC AUC (validatie): {roc_auc_val:.3f}")
    print(f"ROC AUC (test): {roc_auc_test:.3f}")

    # crossvalidation
    k_folds = StratifiedKFold(n_splits=10, shuffle=True)
    scores_cross_validation = cross_val_score(gb, train_x_resampled, train_y_resampled, cv=k_folds)

    print("Random forest tree")
    print("Average CV Score: ", scores_cross_validation.mean())
    print("Number of CV Scores used in Average: ", len(scores_cross_validation))
    print("")

    # feature importance
    result = permutation_importance(gb, val_x_labels, val_y_labels, n_repeats=2, random_state=42)
    sorted_idx = result.importances_mean.argsort()[::-1]
    feature_names = ['ani', "correlation", 'fragmentation', 'average contig length', 'contig_count']
    for i in sorted_idx[:10]:
        print(f"{feature_names[i]}: {result.importances_mean[i]:.4f}")
    return gb_accuracy, gb_f1, mcc_gb


def xgboost_model(train_x_labels, train_y_labels, test_x_labels, test_y_labels, val_x_labels,
                                      val_y_labels):
    """
    Trains and tests  the XGBoost model. Also calculates the accuracy, f1 score
    and MCC score of the model and plots a confusion matrix
    :param train_x_labels: pandas dataframe with the Feature set for training.
    :param train_y_labels: pandas dataframe with the Labels for training data.
    :param test_x_labels: pandas dataframe with the Feature set for testing.
    :param test_y_labels: pandas dataframe with the Labels for testing data.
    :param val_x_labels: pandas dataframe with the rows which will be used to test the model.
    :param val_y_labels: pandas dataframe with the label of contamination for each row in the test_x_label
            dataframe.
    :return:
        accuracy_xgb (float): shows how much contamination the model predicted correctly.
        xgb_f1 (float): a combination of precision and recall which is used for evaluating models.
        mcc_xgb (float): a score between -1 and 1 which is used for evaluating models.
    """
    print("xgboost model")
    labels = [0, 1]

    # Feature selectie
    train_x_labels = train_x_labels.iloc[:, -5:].values.astype(np.float32)
    test_x_labels = test_x_labels.iloc[:, -5:].values.astype(np.float32)
    val_x_labels = val_x_labels.iloc[:, -5:].values.astype(np.float32)

    # Labels
    train_y_labels = train_y_labels.astype(int)
    test_y_labels = test_y_labels.astype(int)
    val_y_labels = val_y_labels.astype(int)

    smote = SMOTE(sampling_strategy=0.1, random_state=42)  # 10% of dataset is minority class
    train_x_resampled, train_y_resampled = smote.fit_resample(train_x_labels, train_y_labels)

    # scale
    scaler = MinMaxScaler()
    train_x_resampled = scaler.fit_transform(train_x_resampled)
    test_x_labels = scaler.transform(test_x_labels)
    val_x_labels = scaler.transform(val_x_labels)

    # train model
    xgb_model = xgb.XGBClassifier(
        n_estimators=10000,
        min_child_weight=2,
        subsample=0.8,
        learning_rate=0.01,
        enable_categorical=True,
        max_delta_step=2,
        objective="binary:logistic",
        gamma=0.1,
        colsample_bytree=0.8,
        random_state=42
    )

    skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
    scores = cross_val_score(xgb_model, train_x_resampled, train_y_resampled, cv=skf, scoring='f1_weighted')
    print("Gemiddelde cross val score:", scores.mean())
    xgb_model.fit(train_x_resampled, train_y_resampled)

    # === VALIDATIE EVALUATIE ===
    val_probs = xgb_model.predict_proba(val_x_labels)[:, 1]
    fpr, tpr, thresholds = roc_curve(val_y_labels, val_probs)
    optimal_idx = np.argmax(tpr - fpr)
    optimal_threshold = thresholds[optimal_idx]
    print(f"Optimale threshold (gebaseerd op validatie): {optimal_threshold:.2f}")

    val_preds = (val_probs > optimal_threshold).astype(int)

    print("=== Validatieset ===")
    print(classification_report(val_y_labels, val_preds, digits=3))
    cm_val = confusion_matrix(val_y_labels, val_preds)
    ConfusionMatrixDisplay(confusion_matrix=cm_val, display_labels=labels).plot()
    plt.title("Confusion Matrix - Validatie - xg")
    plt.show()

    # === TEST EVALUATIE ===
    test_probs = xgb_model.predict_proba(test_x_labels)[:, 1]
    test_preds = (test_probs > optimal_threshold).astype(int)

    print("=== Testset ===")
    print(classification_report(test_y_labels, test_preds, digits=3))
    cm_test = confusion_matrix(test_y_labels, test_preds)
    ConfusionMatrixDisplay(confusion_matrix=cm_test, display_labels=labels).plot()
    plt.title("Confusion Matrix - Test -xg")
    plt.show()

    # VALIDATIE: MCC en F1
    val_f1 = f1_score(val_y_labels, val_preds, average='weighted')
    val_mcc = matthews_corrcoef(val_y_labels, val_preds)
    print(f"Validation F1-score (weighted): {val_f1:.3f}")
    print(f"Validation MCC: {val_mcc:.3f}")

    # TEST: MCC en F1
    test_f1 = f1_score(test_y_labels, test_preds, average='weighted')
    test_mcc = matthews_corrcoef(test_y_labels, test_preds)
    print(f"Test F1-score (weighted): {test_f1:.3f}")
    print(f"Test MCC: {test_mcc:.3f}")

    xgb_f1 = f1_score(test_y_labels, test_preds, average="weighted")
    accuracy_xgb = accuracy_score(test_y_labels, test_preds)
    mcc_xgb = matthews_corrcoef(test_y_labels, test_preds)

    # ROC AUC
    roc_auc_val = roc_auc_score(val_y_labels, val_probs)
    roc_auc_test = roc_auc_score(test_y_labels, test_probs)

    fpr_val, tpr_val, _ = roc_curve(val_y_labels, val_probs)
    fpr_test, tpr_test, _ = roc_curve(test_y_labels, test_probs)

    plt.figure(figsize=(8, 6))
    plt.plot(fpr_val, tpr_val, label=f'Validatie ROC curve (AUC = {roc_auc_val:.2f})')
    plt.plot(fpr_test, tpr_test, label=f'Test ROC curve (AUC = {roc_auc_test:.2f})', linestyle='--')
    plt.plot([0, 1], [0, 1], color='gray', linestyle=':', label='Random Classifier')

    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve - xgb')
    plt.legend(loc='lower right')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    print(f"ROC AUC (validatie): {roc_auc_val:.3f}")
    print(f"ROC AUC (test): {roc_auc_test:.3f}")

    # crossvalidation
    k_folds = StratifiedKFold(n_splits=10, shuffle=True)
    scores_cross_validation = cross_val_score(xgb_model, train_x_resampled, train_y_resampled, cv=k_folds)

    print("Random forest tree")
    print("Average CV Score: ", scores_cross_validation.mean())
    print("Number of CV Scores used in Average: ", len(scores_cross_validation))
    print("")

    # feature importance
    result = permutation_importance(xgb_model, val_x_labels, val_y_labels, n_repeats=10, random_state=42)
    sorted_idx = result.importances_mean.argsort()[::-1]
    feature_names = ['ani', "correlation", 'fragmentation', 'average contig length', 'contig_count']
    for i in sorted_idx[:10]:
        print(f"{feature_names[i]}: {result.importances_mean[i]:.4f}")
    return accuracy_xgb, xgb_f1, mcc_xgb


def train_model_k_nearest_neighbour(train_x_labels, train_y_labels, test_x_labels, test_y_labels, val_x_labels,
                                      val_y_labels):
    """
    Trains and tests K Nearest Neighbours model. Also calculates the accuracy, f1 score
    and MCC score of the model and plots a confusion matrix.
    :param train_x_labels: pandas dataframe with the Feature set for training.
    :param train_y_labels: pandas dataframe with the Labels for training data.
    :param test_x_labels: pandas dataframe with the Feature set for testing.
    :param test_y_labels: pandas dataframe with the Labels for testing data.
    :param val_x_labels: pandas dataframe with the rows which will be used to test the model.
    :param val_y_labels: pandas dataframe with the label of contamination for each row in the test_x_label
            dataframe.
    :return:
        accuracy_score_knn (float): shows how much contamination the model predicted correctly.
        KNN_f1 (float): a combination of precision and recall which is used for evaluating models.
        mcc_knn (float): a score between -1 and 1 which is used for evaluating models.
    """
    labels = [0, 1]
    print("K Nearest Neighbour")

    # Feature selectie
    train_x_labels = train_x_labels.iloc[:, -5:].values.astype(np.float32)
    test_x_labels = test_x_labels.iloc[:, -5:].values.astype(np.float32)
    val_x_labels = val_x_labels.iloc[:, -5:].values.astype(np.float32)

    # Labels
    train_y_labels = train_y_labels.astype(int)
    test_y_labels = test_y_labels.astype(int)
    val_y_labels = val_y_labels.astype(int)

    # oversample minority class
    smote = SMOTE(sampling_strategy=0.15, random_state=42)  # 15% of dataset is minority class
    train_x_resampled, train_y_resampled = smote.fit_resample(train_x_labels, train_y_labels)

    # Scaling
    scaler = MinMaxScaler()
    train_x_resampled = scaler.fit_transform(train_x_resampled)
    test_x_labels = scaler.transform(test_x_labels)
    val_x_labels = scaler.transform(val_x_labels)

    # train model
    knn = KNeighborsClassifier(n_neighbors=100, weights="distance")
    knn.fit(train_x_resampled, train_y_resampled)

    # === VALIDATIE EVALUATIE ===
    val_probs = knn.predict_proba(val_x_labels)[:, 1]
    fpr, tpr, thresholds = roc_curve(val_y_labels, val_probs)
    optimal_idx = np.argmax(tpr - fpr)
    optimal_threshold = thresholds[optimal_idx]
    print(f"Optimale threshold (gebaseerd op validatie): {optimal_threshold:.2f}")

    val_preds = (val_probs > optimal_threshold).astype(int)

    print("=== Validatieset ===")
    print(classification_report(val_y_labels, val_preds, digits=3))
    cm_val = confusion_matrix(val_y_labels, val_preds)
    ConfusionMatrixDisplay(confusion_matrix=cm_val, display_labels=labels).plot()
    plt.title("Confusion Matrix - Validatie - knn")
    plt.show()

    # === TEST EVALUATIE ===
    test_probs = knn.predict_proba(test_x_labels)[:, 1]
    test_preds = (test_probs > optimal_threshold).astype(int)

    print("=== Testset ===")
    print(classification_report(test_y_labels, test_preds, digits=3))
    cm_test = confusion_matrix(test_y_labels, test_preds)
    ConfusionMatrixDisplay(confusion_matrix=cm_test, display_labels=labels).plot()
    plt.title("Confusion Matrix - Test - knn")
    plt.show()

    # VALIDATIE: MCC en F1
    val_f1 = f1_score(val_y_labels, val_preds, average='weighted')
    val_mcc = matthews_corrcoef(val_y_labels, val_preds)
    print(f"Validation F1-score (weighted): {val_f1:.3f}")
    print(f"Validation MCC: {val_mcc:.3f}")

    # TEST: MCC en F1
    test_f1 = f1_score(test_y_labels, test_preds, average='weighted')
    test_mcc = matthews_corrcoef(test_y_labels, test_preds)
    print(f"Test F1-score (weighted): {test_f1:.3f}")
    print(f"Test MCC: {test_mcc:.3f}")

    # Calculates the accuracy, f1 and mcc scores.
    accuracy_score_knn = accuracy_score(test_y_labels, test_preds)
    KNN_f1 = f1_score(test_y_labels, test_preds, average="weighted")
    mcc_knn = matthews_corrcoef(test_y_labels, test_preds)

    # ROC AUC
    roc_auc_val = roc_auc_score(val_y_labels, val_probs)
    roc_auc_test = roc_auc_score(test_y_labels, test_probs)

    fpr_val, tpr_val, _ = roc_curve(val_y_labels, val_probs)
    fpr_test, tpr_test, _ = roc_curve(test_y_labels, test_probs)

    plt.figure(figsize=(8, 6))
    plt.plot(fpr_val, tpr_val, label=f'Validatie ROC curve (AUC = {roc_auc_val:.2f})')
    plt.plot(fpr_test, tpr_test, label=f'Test ROC curve (AUC = {roc_auc_test:.2f})', linestyle='--')
    plt.plot([0, 1], [0, 1], color='gray', linestyle=':', label='Random Classifier')

    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve - knn')
    plt.legend(loc='lower right')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    print(f"ROC AUC (validatie): {roc_auc_val:.3f}")
    print(f"ROC AUC (test): {roc_auc_test:.3f}")

    # crossvalidation
    k_folds = StratifiedKFold(n_splits=10, shuffle=True)
    scores_cross_validation = cross_val_score(knn, train_x_resampled, train_y_resampled, cv=k_folds)

    print("Average CV Score: ", scores_cross_validation.mean())
    print("Number of CV Scores used in Average: ", len(scores_cross_validation))
    print("")
    return accuracy_score_knn, KNN_f1, mcc_knn


def ensemble_classifier_model(train_x_labels, train_y_labels, test_x_labels, test_y_labels, val_x_labels,
                                      val_y_labels):
    """
    Trains and tests the ensemble model. Validation is done with k-fold cross validation.
    Also calculates the accuracy, f1 score  and MCC score of the model and plots a confusion matrix.
    :param train_x_labels: pandas dataframe with the Feature set for training.
    :param train_y_labels: pandas dataframe with the Labels for training data.
    :param test_x_labels: pandas dataframe with the Feature set for testing.
    :param test_y_labels: pandas dataframe with the Labels for testing data.
    :param val_x_labels: pandas dataframe with the rows which will be used to test the model.
    :param val_y_labels: pandas dataframe with the label of contamination for each row in the test_x_label
            dataframe.
    :return:
        accuracy_score_ensemble (float): shows how much contamination the model predicted correctly.
        f1_score_ensemble (float): a combination of precision and recall which is used for evaluating models.
        mcc_vo (float): a score between -1 and 1 which is used for evaluating models.
    """
    print("ensemble classifier")
    labels = [0, 1]

    # Feature selectie
    train_x_labels = train_x_labels.iloc[:, -5:].values.astype(np.float32)
    test_x_labels = test_x_labels.iloc[:, -5:].values.astype(np.float32)
    val_x_labels = val_x_labels.iloc[:, -5:].values.astype(np.float32)

    # Labels
    train_y_labels = train_y_labels.astype(int)
    test_y_labels = test_y_labels.astype(int)
    val_y_labels = val_y_labels.astype(int)

    # oversample minority class
    smote = SMOTE(sampling_strategy=0.15, random_state=42)  # Verhouding 15% van klasse 1
    train_x_resampled, train_y_resampled = smote.fit_resample(train_x_labels, train_y_labels)

    # give weight based on label
    class_weights = compute_class_weight("balanced", classes=np.array([0, 1]), y=train_y_resampled)
    class_weight_dict = {0: class_weights[0], 1: class_weights[1]}

    # Scaling
    scaler = MinMaxScaler()
    train_x_resampled = scaler.fit_transform(train_x_resampled)
    test_x_labels = scaler.transform(test_x_labels)
    val_x_labels = scaler.transform(val_x_labels)

    # train models
    rf = RandomForestClassifier(n_estimators=10000,
                                max_features='sqrt',
                                class_weight=class_weight_dict,
                                bootstrap=True,
                                min_samples_split=5,
                                min_samples_leaf=5,
                                random_state=42)
    et = ExtraTreesClassifier(n_estimators=10000,
                              max_features='sqrt',
                              class_weight=class_weight_dict)
    xg = xgb.XGBClassifier(n_estimators=10000,
                            min_child_weight=2,
                            subsample=0.8,
                            learning_rate=0.01,
                            enable_categorical=True,
                            max_delta_step=2,
                            objective="binary:logistic",
                            gamma=0.1,
                            colsample_bytree=0.8,
                            random_state=42)
    gb = GradientBoostingClassifier(n_estimators=10000,
                                    learning_rate=0.1,
                                    subsample=0.8,
                                    max_features=None)
    ba = BaggingClassifier(n_estimators=10000)

    ensemble = VotingClassifier(estimators=[("rf", rf), ("et", et), ("xg", xg), ("gb", gb), ("ba", ba)], voting="soft",)

    k_folds = StratifiedKFold(n_splits=3, shuffle=True)
    scores = cross_val_score(ensemble, train_x_resampled, train_y_resampled, cv=k_folds)
    print("Cross Validation Scores: ", scores)
    print("Average CV Score: ", scores.mean())
    print("Number of CV Scores used in Average: ", len(scores))
    print("")
    # Plots the score of every k-fold.
    plt.plot(range(len(scores)), scores)
    plt.xlabel("K-folds")
    plt.ylabel("Cross validation score")
    plt.title("Cross Validation scores per k-fold")
    plt.show()

    # Fits the model onto the training data.
    ensemble.fit(train_x_resampled, train_y_resampled)

    # === VALIDATIE EVALUATIE ===
    val_probs = ensemble.predict_proba(val_x_labels)[:, 1]
    fpr, tpr, thresholds = roc_curve(val_y_labels, val_probs)
    optimal_idx = np.argmax(tpr - fpr)
    optimal_threshold = thresholds[optimal_idx]
    print(f"Optimale threshold (gebaseerd op validatie): {optimal_threshold:.2f}")

    val_preds = (val_probs > optimal_threshold).astype(int)

    print("=== Validatieset ===")
    print(classification_report(val_y_labels, val_preds, digits=3))
    cm_val = confusion_matrix(val_y_labels, val_preds)
    ConfusionMatrixDisplay(confusion_matrix=cm_val, display_labels=labels).plot()
    plt.title("Confusion Matrix - Validatie - ensemble")
    plt.show()

    # === TEST EVALUATIE ===
    test_probs = ensemble.predict_proba(test_x_labels)[:, 1]
    test_preds = (test_probs > optimal_threshold).astype(int)

    print("=== Testset ===")
    print(classification_report(test_y_labels, test_preds, digits=3))
    cm_test = confusion_matrix(test_y_labels, test_preds)
    ConfusionMatrixDisplay(confusion_matrix=cm_test, display_labels=labels).plot()
    plt.title("Confusion Matrix - Test - ensemble")
    plt.show()

    # VALIDATIE: MCC en F1
    val_f1 = f1_score(val_y_labels, val_preds, average='weighted')
    val_mcc = matthews_corrcoef(val_y_labels, val_preds)
    print(f"Validation F1-score (weighted): {val_f1:.3f}")
    print(f"Validation MCC: {val_mcc:.3f}")

    # TEST: MCC en F1
    test_f1 = f1_score(test_y_labels, test_preds, average='weighted')
    test_mcc = matthews_corrcoef(test_y_labels, test_preds)
    print(f"Test F1-score (weighted): {test_f1:.3f}")
    print(f"Test MCC: {test_mcc:.3f}")

    # Calculates the accuracy, f1 and mcc score.
    accuracy_score_ensemble = accuracy_score(test_y_labels, test_preds)
    f1_score_ensemble = f1_score(test_y_labels, test_preds, average="weighted")
    mcc_vo = matthews_corrcoef(test_y_labels, test_preds)
    print("Accuracy score: ", (accuracy_score_ensemble / 1) * 100, "%")
    print("F1 score: ", f1_score_ensemble)
    print("MCC: ", mcc_vo)
    print("precision: ", precision_score(test_y_labels, test_preds))
    print("")

    # ROC AUC
    roc_auc_val = roc_auc_score(val_y_labels, val_probs)
    roc_auc_test = roc_auc_score(test_y_labels, test_probs)

    fpr_val, tpr_val, _ = roc_curve(val_y_labels, val_probs)
    fpr_test, tpr_test, _ = roc_curve(test_y_labels, test_probs)

    plt.figure(figsize=(8, 6))
    plt.plot(fpr_val, tpr_val, label=f'Validatie ROC curve (AUC = {roc_auc_val:.2f})')
    plt.plot(fpr_test, tpr_test, label=f'Test ROC curve (AUC = {roc_auc_test:.2f})', linestyle='--')
    plt.plot([0, 1], [0, 1], color='gray', linestyle=':', label='Random Classifier')

    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve - ensemble')
    plt.legend(loc='lower right')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    print(f"ROC AUC (validatie): {roc_auc_val:.3f}")
    print(f"ROC AUC (test): {roc_auc_test:.3f}")

    # crossvalidation
    k_folds = StratifiedKFold(n_splits=10, shuffle=True)
    scores_cross_validation = cross_val_score(ensemble, train_x_resampled, train_y_resampled, cv=k_folds)

    print("ensemble")
    print("Average CV Score: ", scores_cross_validation.mean())
    print("Number of CV Scores used in Average: ", len(scores_cross_validation))
    print("")

    # feature importance
    result = permutation_importance(ensemble, val_x_labels, val_y_labels, n_repeats=10, random_state=42)
    sorted_idx = result.importances_mean.argsort()[::-1]
    feature_names = ['ani', "correlation", 'fragmentation', 'average contig length', 'contig_count']
    for i in sorted_idx[:10]:
        print(f"{feature_names[i]}: {result.importances_mean[i]:.4f}")

    # Saves the model into a file.
    joblib.dump(ensemble, "ensemble_2xaugment_mcc0.54.joblib")
    return accuracy_score_ensemble, f1_score_ensemble, mcc_vo


def plot_accuracy(accuracy_score_random_forest, accuracy_score_KNN,
                  accuracy_score_ensemble, gb_accuracy, accuracy_xgb, accuracy_score_mlp,
                  accuracy_score_Naive_bayes, accuracy_score_complementNaive_bayes):
    """
    Makes the figure showing the f1 scores of the different models.

    Args:
        accuracy_score_random_forest (float): shows how much contamination the random forest model predicted correctly.
        accuracy_score_KNN (float): shows how much contamination the KNN model predicted correctly.
        accuracy_score_Naive_bayes (float): shows how much contamination the Naive Bayes model predicted correctly.
        accuracy_score_lda (float): shows how much contamination the lda model predicted correctly.
        accuracy_score_qda (float): shows how much contamination the qda model predicted correctly.
        accuracy_score_ensemble (float): shows how much contamination the ensemble model predicted correctly.
        accuracy_xgb (float): shows how much contamination the XGBoost model predicted correctly.
        accuracy_score_mlp(float): shows how much contamination the mlp model predicted correctly.
        gb_accuracy(float): shows how much contamination the gradient boost model predicted correctly.
        accuracy_score_complementNaive_bayes(float): shows how much contamination the complement naive bayes model predicted correctly.
    """
    # Gets the accuracy scores from the models and puts them into a list.
    list_score = [accuracy_score_random_forest, accuracy_score_KNN,  accuracy_score_ensemble, gb_accuracy, accuracy_xgb,
                  accuracy_score_mlp, accuracy_score_Naive_bayes, accuracy_score_complementNaive_bayes]
    models = ["Random Forest", "K nearest neighbor", "Voting classifier", "Gradient boost", "XG boost",
              "multilayer perceptron", "accuracy_score Naive_bayes", "complementNaive_bayes"]

    # Creating the bar plot.
    plt.bar(models, list_score, color="blue", width=0.4)
    plt.xticks(rotation=12)
    plt.xlabel("Models")
    plt.ylabel("Accuracy")
    plt.title("Accuracy score of the models")
    plt.show()


def plot_f1(randomforest_tree_f1, KNN_f1,  f1_score_ensemble, gb_f1, xgb_f1, mlp_f1, f1_qda, f1_lda, Naive_bayes_f1,
            complementNaive_bayes_f1):
    """
    Makes the figure showing the f1 scores of the different models.
    Args:
        randomforest_tree_f1 (float): a combination of precision and recall which is used for evaluating the rf model.
        KNN_f1 (float): a combination of precision and recall which is used for evaluating the knn model.
        Naive_bayes_f1 (float): a combination of precision and recall which is used for evaluating the Naive bayes
        model.
        f1_lda (float): a combination of precision and recall which is used for evaluating the lda model.
        f1_qda (float): a combination of precision and recall which is used for evaluating the qda model.
        f1_score_ensemble (float): a combination of precision and recall which is used for evaluating the ensemble
        model.
        gb_f1 (float): a combination of precision and recall which is used for evaluating the gb model.
        xgb_f1 (float): a combination of precision and recall which is used for evaluating the xgb model.
        mlp_f1 (float): a combination of precision and recall which is used for evaluating the xgb model.
        complementNaive_bayes_f1 (float): a combination of precision and recall which is used for evaluating the xgb
        model.
    """
    # Gets the accuracy scores from the models and puts them into a list.
    list_score = [randomforest_tree_f1, KNN_f1, f1_score_ensemble, gb_f1, xgb_f1, mlp_f1,  f1_qda, f1_lda,
                  Naive_bayes_f1, complementNaive_bayes_f1]
    models = ["Random Forest", "K nearest neighbor",
              "Voting classifier", "Gradient boost", "XG boost", "multilayer perceptron", "qda", "lda", "naive",
              "complement"]

    # Creating the bar plot.
    plt.bar(models, list_score, color="blue", width=0.4)
    plt.xticks(rotation=12)
    plt.xlabel("Models")
    plt.ylabel("F1 score")
    plt.title("F1 score of the models")
    plt.show()


def plot_mcc(mcc_rf, mcc_knn, mcc_vo, mcc_gb, mcc_xgb, mcc_mlp, mcc_qda, mcc_lda, mcc_naive, mcc_complementnaive):
    """
    Makes the figure showing the mcc scores of the different models.
    Args:
        mcc_rf (float): a score between -1 and 1 which is used for evaluating the random forest model.
        mcc_knn (float): a score between -1 and 1 which is used for evaluating the knn model.
        mcc_naive (float): a score between -1 and 1 which is used for evaluating the naive bayes model.
        mcc_lda (float): a score between -1 and 1 which is used for evaluating the lda model.
        mcc_qda (float): a score between -1 and 1 which is used for evaluating the qda model.
        mcc_vo (float): a score between -1 and 1 which is used for evaluating the ensemble model.
        mcc_gb (float): a score between -1 and 1 which is used for evaluating the gb model.
        mcc_xgb (float): a score between -1 and 1 which is used for evaluating the xgb model.
        mcc_mlp (float): a score between -1 and 1 which is used for evaluating the xgb model.
        mcc_complementnaive (float): a score between -1 and 1 which is used for evaluating the xgb model.
    """
    # Gets the accuracy scores from the models and puts them into a list.
    list_score = [mcc_rf, mcc_knn, mcc_vo, mcc_gb,  mcc_xgb, mcc_mlp, mcc_qda, mcc_lda, mcc_naive, mcc_complementnaive]
    models = ["Random Forest", "K nearest neighbor",
              "Voting classifier", "Gradient boost", "XG boost", "multilayer perceptron", "qda", "lda", "naive",
              "complement"]

    # Creating the bar plot.
    plt.bar(models, list_score, color="blue", width=0.4)
    plt.xticks(rotation=12)
    plt.xlabel("Models")
    plt.ylabel("MCC score")
    plt.title("MCC score of the models")
    plt.show()


if __name__ == '__main__':
    # input
    train_file = "..\\..\\output_from_scripts\\test_train_data_ml\\train_set.csv"
    test_file = "..\\..\\output_from_scripts\\test_train_data_ml\\test_set.csv"
    val_file = "..\\..\\output_from_scripts\\test_train_data_ml\\val_set.csv"

    test_data, train_data, val_data = read_file(train_file, test_file, val_file)

    train_x_labels, train_y_labels, test_x_labels, test_y_labels, val_x_labels, val_y_labels \
        = get_labels(train_data, test_data, val_data)

    accuracy_score_complementNaive_bayes, complementNaive_bayes_f1, mcc_complementnaive = \
        train_model_complementnaive_bayes(train_x_labels, train_y_labels, test_x_labels, test_y_labels, val_x_labels,
                                          val_y_labels)

    accuracy_score_random_forest, randomforest_tree_f1, mcc_rf = (
        train_model_random_forest(train_x_labels, train_y_labels, test_x_labels, test_y_labels,
                                  val_x_labels, val_y_labels))
    gb_accuracy, gb_f1, mcc_gb = gradient_boost(train_x_labels, train_y_labels, test_x_labels, test_y_labels,
                                                val_x_labels, val_y_labels)

    accuracy_xgb, xgb_f1, mcc_xgb = xgboost_model(train_x_labels, train_y_labels, test_x_labels, test_y_labels,
                                                  val_x_labels, val_y_labels)

    accuracy_score_KNN, KNN_f1, mcc_knn = train_model_k_nearest_neighbour(train_x_labels, train_y_labels,
                                                                          test_x_labels, test_y_labels, val_x_labels,
                                                                          val_y_labels)
    accuracy_score_ensemble, f1_score_ensemble, mcc_vo = ensemble_classifier_model(train_x_labels, train_y_labels,
                                                                                   test_x_labels, test_y_labels,
                                                                                   val_x_labels, val_y_labels)
    accuracy_score_mlp, mlp_f1, mcc_mlp = train_model_mlp(train_x_labels, train_y_labels, test_x_labels, test_y_labels,
                                                          val_x_labels, val_y_labels)
    accuracy_score_Naive_bayes, Naive_bayes_f1, mcc_naive = train_model_naive_bayes(train_x_labels, train_y_labels,
                                                                                    test_x_labels, test_y_labels,
                                                                                    val_x_labels, val_y_labels)

    accuracy_score_lda, f1_lda, mcc_lda = linear_discriminant_analysis_model(train_x_labels, train_y_labels,
                                                                             test_x_labels, test_y_labels, val_x_labels,
                                                                             val_y_labels)

    accuracy_score_qda, f1_qda, mcc_qda = quadratic_discriminant_model(train_x_labels, train_y_labels, test_x_labels,
                                                                       test_y_labels, val_x_labels,
                                                                       val_y_labels)
    plot_accuracy(accuracy_score_random_forest, accuracy_score_KNN, accuracy_score_ensemble, gb_accuracy, accuracy_xgb,
                  accuracy_score_mlp, accuracy_score_Naive_bayes, accuracy_score_complementNaive_bayes)
    plot_f1(randomforest_tree_f1, KNN_f1,  f1_score_ensemble, gb_f1, xgb_f1, mlp_f1, f1_qda, f1_lda, Naive_bayes_f1,
            complementNaive_bayes_f1)
    plot_mcc(mcc_rf, mcc_knn, mcc_vo, mcc_gb, mcc_xgb, mcc_mlp, mcc_qda, mcc_lda, mcc_naive, mcc_complementnaive)