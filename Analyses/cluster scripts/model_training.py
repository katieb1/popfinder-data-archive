import pandas as pd
from popfinder.classifier import PopClassifier
import argparse

parser = argparse.ArgumentParser(description='PopFinder Training')
parser.add_argument('-o', type=str, help='Output folder from PopFinder Tuning Script')
args = vars(parser.parse_args())
output_folder = args['o']

# Testing 
output_folder = "F:/PopFinder"

# Load hyperparameter tuning results and remove performance metrics
results = pd.read_csv(f"{output_folder}/hyperparam_search_results.csv")
best_hyperparams = results.loc[results['valid_mcc'].idxmax()]
hp_subset = len(best_hyperparams) - 7 # 7 is number of performance metrics
best_hyperparams = best_hyperparams[0:hp_subset]
best_hyperparams = best_hyperparams.to_dict()

# Change number of epochs to 2000
best_hyperparams['epochs'] = 2000

# Load classifier object
classifier = PopClassifier.load(f"{output_folder}/classifier.pkl")
classifier.output_folder = f"{output_folder}/popfinder_train"

# Train using best hyperparameters
classifier.train(cv_splits=5, nreps=5, bootstraps=5, jobs=1,
                 **best_hyperparams)

classifier.save()