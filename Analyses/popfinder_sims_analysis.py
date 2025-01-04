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
sim_dir = "F:/PopFinder/Results/popfinder/"
LS_folder = sim_dir + "LS_outputs/"
LS_train = LS_folder + "popfinder_train/"
MS_folder = sim_dir + "MS_outputs/"
MS_train = MS_folder + "popfinder_train/"
HS_folder = sim_dir + "HS_outputs/"
HS_train = HS_folder + "popfinder_train/"

#%% Function for summarizing assignment performance metrics
def get_performance_metrics(cm, y_true, y_pred, labels, runtime, digits=3):

    accuracy = np.round(cm.diagonal(), digits)
    precision = np.round(np.mean(cm.diagonal()/cm.sum(axis=0)), digits)
    recall = np.round(np.mean(cm.diagonal()/cm.sum(axis=1)), digits)
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

#%% Weak Structure ------------------------------------------------------------
#%% Load classifier
LSclassifier = PopClassifier.load(os.path.join(LS_folder, "classifier.pkl"))
LSclassifier.output_folder = LS_folder

#%% Test classifier - using best model
LSclassifier.test(use_best_model=True)
LSclassifier.plot_confusion_matrix()
LSclassifier.get_test_summary()

#%% Test classifier - using all models 
LSclassifier.test(use_best_model=False, ensemble_accuracy_threshold=0.5)
LSclassifier.plot_confusion_matrix()
LSclassifier.get_test_summary()

#%% Get most important SNPs
LSclassifier.rank_site_importance()

#%% Use classifier to assign populations to unknowns
start = time.time()
unknown_assignments = LSclassifier.assign_unknown(use_best_model=False, ensemble_accuracy_threshold=0.5)
end = time.time()
runtime = end - start

x = np.array(["p1", "p2", "p3", "p4", "p5"])
real_pops = np.repeat(x, 960)
unknown_assignments["real_pop"] = real_pops
unknown_assignments["success"] = unknown_assignments["most_assigned_pop_across_models"] == unknown_assignments["real_pop"]

overall_acc = sum(unknown_assignments["success"]) / len(unknown_assignments)

# %%
# Create confusion matrix
conf_matrix = confusion_matrix(unknown_assignments["real_pop"], 
                               unknown_assignments["most_assigned_pop_across_models"],
                               labels = x, normalize="true")
disp = ConfusionMatrixDisplay(confusion_matrix=conf_matrix,
                              display_labels=x)
disp.plot()
plt.show()

#%% Get performance metrics for assignment
perf_metrics = get_performance_metrics(
    conf_matrix, unknown_assignments["real_pop"],
    unknown_assignments["most_assigned_pop_across_models"], x,
    runtime)

perf_metrics.to_csv(os.path.join(LS_folder, "unknown_metrics.csv"), index=False)


# %% Intermediate Structure -----------------------------------------------------
#%% Load classifier
MSclassifier = PopClassifier.load(os.path.join(MS_folder, "classifier.pkl"))
MSclassifier.output_folder = MS_folder

#%% Test classifier - using best model
MSclassifier.test(use_best_model=True)
MSclassifier.plot_confusion_matrix()
MSclassifier.get_test_summary()

#%% Test classifier - using all models 
MSclassifier.test(use_best_model=False)
MSclassifier.plot_confusion_matrix()
MSclassifier.get_test_summary()

#%% Get most important SNPs
MSclassifier.rank_site_importance()

#%% Use classifier to assign populations to unknowns
start = time.time()
unknown_assignments = MSclassifier.assign_unknown(use_best_model=False)
end = time.time()
runtime = end - start

x = np.array(["p1", "p2", "p3", "p4", "p5"])
real_pops = np.repeat(x, 960)
unknown_assignments["real_pop"] = real_pops
unknown_assignments["success"] = unknown_assignments["most_assigned_pop_across_models"] == unknown_assignments["real_pop"]
overall_acc = sum(unknown_assignments["success"]) / len(unknown_assignments)

# %%
# Create confusion matrix
conf_matrix = confusion_matrix(unknown_assignments["real_pop"], 
                               unknown_assignments["most_assigned_pop_across_models"],
                               labels = x, normalize="true")
disp = ConfusionMatrixDisplay(confusion_matrix=conf_matrix,
                              display_labels=x)
disp.plot()
plt.show()

#%% Get performance metrics for assignment
perf_metrics = get_performance_metrics(
    conf_matrix, unknown_assignments["real_pop"],
    unknown_assignments["most_assigned_pop_across_models"], x,
    runtime)

perf_metrics.to_csv(os.path.join(MS_folder, "unknown_metrics.csv"), index=False)

# %% Strong Structure -----------------------------------------------------
#%% Load classifier
HSclassifier = PopClassifier.load(os.path.join(HS_folder, "classifier.pkl"))
HSclassifier.output_folder = HS_folder

#%% Test classifier - using best model
HSclassifier.test(use_best_model=True)
HSclassifier.plot_confusion_matrix()
HSclassifier.get_test_summary()

#%% Test classifier - using all models 
HSclassifier.test(use_best_model=False)
HSclassifier.plot_confusion_matrix()
HSclassifier.get_test_summary()

#%% Get most important SNPs
HSclassifier.rank_site_importance()

#%% Use classifier to assign populations to unknowns
start = time.time()
unknown_assignments = HSclassifier.assign_unknown(use_best_model=False)
end = time.time()
runtime = end - start

x = np.array(["p1", "p2", "p3", "p4", "p5"])
real_pops = np.repeat(x, 960)
unknown_assignments["real_pop"] = real_pops
unknown_assignments["success"] = unknown_assignments["most_assigned_pop_across_models"] == unknown_assignments["real_pop"]
overall_acc = sum(unknown_assignments["success"]) / len(unknown_assignments)

# %%
# Create confusion matrix
conf_matrix = confusion_matrix(unknown_assignments["real_pop"], 
                               unknown_assignments["most_assigned_pop_across_models"],
                               labels = x, normalize="true")
disp = ConfusionMatrixDisplay(confusion_matrix=conf_matrix,
                              display_labels=x)
disp.plot()
plt.show()

#%% Get performance metrics for assignment
perf_metrics = get_performance_metrics(
    conf_matrix, unknown_assignments["real_pop"],
    unknown_assignments["most_assigned_pop_across_models"], x,
    runtime)

perf_metrics.to_csv(os.path.join(HS_folder, "unknown_metrics.csv"), index=False)

# %%
