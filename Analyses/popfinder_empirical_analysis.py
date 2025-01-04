## Script for processing results from simulated dataset tuning and training

#%% Imports
from popfinder.classifier import PopClassifier
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, matthews_corrcoef
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import time


#%% Constants
input_dir = "F:/PopFinder/Data/empirical data/"
results_dir = "F:/PopFinder/Results/popfinder/"
LESP_folder = "oct27_lesp/"
TBMU_folder = "oct6_tbmu/"
NOFU6_folder = "nofu6/"
NOFU4_folder = "sep11_nofu/"

#%% Function for summarizing assignment performance metrics
def get_performance_metrics(cm, y_true, y_pred, labels, runtime, digits=3):

    # Accuracy per population
    accuracy = np.round(cm.diagonal(), digits)

    # precision per population
    precision = cm.diagonal()/cm.sum(axis=0)
    precision = np.nan_to_num(precision, 0)

    # recall per population
    recall = cm.diagonal()/cm.sum(axis=1)
    recall = np.nan_to_num(recall, 0)    

    # Overall metrics
    precision = np.round(np.mean(precision), digits)
    recall = np.round(np.mean(recall), digits)
    f1 = np. round(np.mean(2 * (precision * recall) / (precision + recall)), digits)
    mcc = np.round(matthews_corrcoef(y_true, y_pred), digits)

    results_df = pd.DataFrame({
        "Metric": ["Accuracy", "Precision", "Recall", "F1", "MCC"],
        "Value": [np.mean(accuracy), precision, recall, f1, mcc]})
    
    pop_accuracies = pd.DataFrame({
        "Metric": labels,
        "Value": accuracy})
    
    runtime_df = pd.DataFrame({
        "Metric": ["Runtime"],
        "Value": np.round(runtime, digits)})  
    
    results_df = pd.concat([results_df, pop_accuracies, runtime_df], axis=0)

    return results_df

#%% LESP ------------------------------------------------------------
#%% Load classifier
input_folder = os.path.join(input_dir, LESP_folder)
output_folder = os.path.join(results_dir, LESP_folder)
LESPclassifier = PopClassifier.load(os.path.join(output_folder, "classifier.pkl"))
LESPclassifier.output_folder = output_folder

#%% Test classifier - using best model
LESPclassifier.test(use_best_model=True)
LESPclassifier.plot_confusion_matrix()
LESPclassifier.get_test_summary()

#%% Test classifier - using all models 
LESPclassifier.test(use_best_model=False, ensemble_accuracy_threshold=0.25)
LESPclassifier.plot_confusion_matrix()
LESPclassifier.get_test_summary()

#%% Get most influential SNPs
LESPclassifier.rank_site_importance()

#%% Use classifier to assign populations to unknowns
start = time.time()
unknown_assignments = LESPclassifier.assign_unknown(use_best_model=False)
end = time.time()
runtime = end - start

#%% Add real populations 
popmap = pd.read_csv(os.path.join(input_folder, "popmap.txt"), sep="\t")
unknown_assignments = pd.merge(popmap[["sampleID", "pop"]], unknown_assignments, on = "sampleID")
unknown_assignments = unknown_assignments.rename(columns={"pop_x": "real_pop"})
unknown_assignments["success"] = unknown_assignments["most_assigned_pop_across_models"] == unknown_assignments["real_pop"]

overall_acc = sum(unknown_assignments["success"]) / len(unknown_assignments)

# %%
# Create confusion matrix
conf_matrix = confusion_matrix(unknown_assignments["real_pop"], 
                               unknown_assignments["most_assigned_pop_across_models"],
                               labels = np.unique(popmap["pop"]), normalize="true")
disp = ConfusionMatrixDisplay(confusion_matrix=conf_matrix,
                              display_labels=np.unique(popmap["pop"]))
disp.plot()
plt.show()

#%% Get performance metrics for assignment
perf_metrics = get_performance_metrics(
    conf_matrix, unknown_assignments["real_pop"],
    unknown_assignments["most_assigned_pop_across_models"], 
    np.unique(popmap["pop"]),
    runtime)

perf_metrics.to_csv(os.path.join(output_folder, "unknown_metrics.csv"), index=False)

#%% TBMU ------------------------------------------------------------
#%% Load classifier
input_folder = os.path.join(input_dir, TBMU_folder)
output_folder = os.path.join(results_dir, TBMU_folder)
TBMUclassifier = PopClassifier.load(os.path.join(output_folder, "classifier.pkl"))
TBMUclassifier.output_folder = output_folder

#%% Test classifier - using best model
TBMUclassifier.test(use_best_model=True)
TBMUclassifier.plot_confusion_matrix()
TBMUclassifier.get_test_summary()

#%% Test classifier - using all models 
TBMUclassifier.test(use_best_model=False, ensemble_accuracy_threshold=0.1)
TBMUclassifier.plot_confusion_matrix()
TBMUclassifier.get_test_summary()

#%% Get most influential SNPs
TBMUclassifier.rank_site_importance()

#%% Use classifier to assign populations to unknowns
start = time.time()
unknown_assignments = TBMUclassifier.assign_unknown(use_best_model=False)
end = time.time()
runtime = end - start

#%% Add real populations 
popmap = pd.read_csv(os.path.join(input_folder, "popmap.txt"), sep="\t")
unknown_assignments = pd.merge(popmap[["sampleID", "pop"]], unknown_assignments, on = "sampleID")
unknown_assignments = unknown_assignments.rename(columns={"pop_x": "real_pop"})
unknown_assignments["success"] = unknown_assignments["most_assigned_pop_across_models"] == unknown_assignments["real_pop"]

# Remove samples that were not included in training
unknown_assignments = unknown_assignments[~unknown_assignments['real_pop'].isin(['coa', 'cog', 'fun'])]
poplabels = np.unique(unknown_assignments["real_pop"])

overall_acc = sum(unknown_assignments["success"]) / len(unknown_assignments)

# %%
# Create confusion matrix
conf_matrix = confusion_matrix(unknown_assignments["real_pop"], 
                               unknown_assignments["most_assigned_pop_across_models"],
                               labels = poplabels, normalize="true")
disp = ConfusionMatrixDisplay(confusion_matrix=conf_matrix,
                              display_labels=poplabels)
disp.plot()
plt.show()

#%% Get performance metrics for assignment
perf_metrics = get_performance_metrics(
    conf_matrix, unknown_assignments["real_pop"],
    unknown_assignments["most_assigned_pop_across_models"], 
    poplabels,
    runtime)

perf_metrics.to_csv(os.path.join(output_folder, "unknown_metrics.csv"), index=False)

# %% NOFU (6 Colonies) ------------------------------------------------------------
#%% Load classifier
input_folder = os.path.join(input_dir, "nofu")
output_folder = os.path.join(results_dir, NOFU6_folder)
popmap_filename = "popmap_6cols.txt"
NOFU6classifier = PopClassifier.load(os.path.join(output_folder, "classifier.pkl"))
NOFU6classifier.output_folder = output_folder

#%% Test classifier - using best model
NOFU6classifier.test(use_best_model=True)
NOFU6classifier.plot_confusion_matrix()
NOFU6classifier.get_test_summary()

#%% Test classifier - using all models 
NOFU6classifier.test(use_best_model=False)
NOFU6classifier.plot_confusion_matrix()
NOFU6classifier.get_test_summary()

#%% Use classifier to assign populations to unknowns
start = time.time()
unknown_assignments = NOFU6classifier.assign_unknown(use_best_model=False)
end = time.time()
runtime = end - start

#%% Add real populations 
popmap = pd.read_csv(os.path.join(input_folder, popmap_filename), sep="\t")
unknown_assignments = pd.merge(popmap[["sampleID", "pop"]], unknown_assignments, on = "sampleID")
unknown_assignments = unknown_assignments.rename(columns={"pop_x": "real_pop"})
unknown_assignments["success"] = unknown_assignments["most_assigned_pop_across_models"] == unknown_assignments["real_pop"]

# Remove samples that were not included in training
unknown_assignments = unknown_assignments[~unknown_assignments['real_pop'].isin(['coa', 'cog', 'fun'])]
poplabels = np.unique(unknown_assignments["real_pop"])

overall_acc = sum(unknown_assignments["success"]) / len(unknown_assignments)

# %%
# Create confusion matrix
conf_matrix = confusion_matrix(unknown_assignments["real_pop"], 
                               unknown_assignments["most_assigned_pop_across_models"],
                               labels = poplabels, normalize="true")
disp = ConfusionMatrixDisplay(confusion_matrix=conf_matrix,
                              display_labels=poplabels)
disp.plot()
plt.show()

#%% Get performance metrics for assignment
perf_metrics = get_performance_metrics(
    conf_matrix, unknown_assignments["real_pop"],
    unknown_assignments["most_assigned_pop_across_models"], 
    poplabels,
    runtime)

perf_metrics.to_csv(os.path.join(output_folder, "unknown_metrics.csv"), index=False)

# %% NOFU (4 Colonies) ------------------------------------------------------------
#%% Load classifier
input_folder = os.path.join(input_dir, "nofu")
output_folder = os.path.join(results_dir, NOFU4_folder)
popmap_filename = "popmap_4cols.txt"
NOFU4classifier = PopClassifier.load(os.path.join(output_folder, "classifier.pkl"))
NOFU4classifier.output_folder = output_folder

#%% Test classifier - using best model
NOFU4classifier.test(use_best_model=True)
NOFU4classifier.plot_confusion_matrix()
NOFU4classifier.get_test_summary()

#%% Test classifier - using all models 
NOFU4classifier.test(use_best_model=False, ensemble_accuracy_threshold=0.5)
NOFU4classifier.plot_confusion_matrix()
NOFU4classifier.get_test_summary()

#%% Get most influential SNPs
NOFU4classifier.rank_site_importance()

#%% Use classifier to assign populations to unknowns
start = time.time()
unknown_assignments = NOFU4classifier.assign_unknown(use_best_model=False)
end = time.time()
runtime = end - start

#%% Add real populations 
popmap = pd.read_csv(os.path.join(input_folder, popmap_filename), sep="\t")
unknown_assignments = pd.merge(popmap[["sampleID", "pop"]], unknown_assignments, on = "sampleID")
unknown_assignments = unknown_assignments.rename(columns={"pop_x": "real_pop"})
unknown_assignments["success"] = unknown_assignments["most_assigned_pop_across_models"] == unknown_assignments["real_pop"]

# Remove samples that were not included in training
unknown_assignments = unknown_assignments[~unknown_assignments['real_pop'].isin(['coa', 'cog', 'fun'])]
poplabels = np.unique(unknown_assignments["real_pop"])

overall_acc = sum(unknown_assignments["success"]) / len(unknown_assignments)

# %%
# Create confusion matrix
conf_matrix = confusion_matrix(unknown_assignments["real_pop"], 
                               unknown_assignments["most_assigned_pop_across_models"],
                               labels = poplabels, normalize="true")
disp = ConfusionMatrixDisplay(confusion_matrix=conf_matrix,
                              display_labels=poplabels)
disp.plot()
plt.show()

#%% Get performance metrics for assignment
perf_metrics = get_performance_metrics(
    conf_matrix, unknown_assignments["real_pop"],
    unknown_assignments["most_assigned_pop_across_models"], 
    poplabels,
    runtime)

perf_metrics.to_csv(os.path.join(output_folder, "unknown_metrics.csv"), index=False)
