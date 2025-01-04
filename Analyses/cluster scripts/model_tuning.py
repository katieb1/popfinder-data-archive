from popfinder.dataloader import GeneticData
from popfinder.classifier import PopClassifier
from popfinder.tuning import hyperparam_search
import argparse

parser = argparse.ArgumentParser(description='PopFinder Tuning')
parser.add_argument('-g', type=str, help='Path to genetic data file')
parser.add_argument('-s', type=str, help='Path to sample data file')
parser.add_argument('-o', type=str, help='Output folder')
args = vars(parser.parse_args())

genetic_data = args['g']
sample_data = args['s']
output_folder = args['o']

# Testing
# genetic_data = "C:/projects/A280/simulated data/katie_new/model1_LS/model1_LS_500sample.recode.vcf"
# sample_data = "C:/projects/A280/simulated data/katie_new/popmap_w_NAs.txt"
# output_folder = "C:/temp/popfinder_test"

# Load data object
data_obj = GeneticData(genetic_data=genetic_data, sample_data=sample_data, seed=42)

# Create classifier object
classifier = PopClassifier(data_obj, output_folder=output_folder)

# Save classifier object to output folder
classifier.save(save_path=output_folder, filename="classifier.pkl")

# Hyperparameter options
trials = 25 # for testing
nreps = 5 # for testing
bootstraps = 5 # for testing
epoch_options = [1000] # for testing
learning_options = [0.01, 0.001]
dropout_options = [0.025, 0]
batch_options = [16, 32]
hidden_size_options = [16, 32]
hidden_layers_options = [2, 5]

hyperparam_dict = dict({"optimizer": ["Adam"],
                        "beta1": [0.9, 0.95, 0.99],
                        "beta2": [0.999, 0.9999],
                        "weight_decay": [0, 0.01, 0.1],
                        "epsilon": [1e-8, 1e-7, 1e-6]})

results = hyperparam_search(
    classifier, trials=trials, cv_splits=5, nreps=nreps, 
    bootstraps=bootstraps, patience=10, 
    min_delta=0.001, learning_rate=learning_options, 
    dropout_prop=dropout_options, batch_size=batch_options, 
    hidden_size=hidden_size_options, hidden_layers=hidden_layers_options, 
    epochs=epoch_options, jobs=5, hyperparam_dict=hyperparam_dict)

# Save results
results.to_csv(f"{output_folder}/hyperparam_search_results.csv", index=False)