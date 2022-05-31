#!/usr/bin/env python3

####################
# Import libraries #
####################

import os 
import sys
import datetime
import numpy
import pandas
import scanpy 
import anndata
import scipy 
import matplotlib.pyplot
import functools
import random
import itertools
import warnings
import argparse
import seaborn

#######################
# Import ML libraries #
#######################

# Required libraries regardless of the model you choose

from sklearn.metrics import confusion_matrix, classification_report, accuracy_score, classification_report
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline

# Library for Logistic Regression
from sklearn.linear_model import LogisticRegression

# Library for Random Forest 
from sklearn.ensemble import RandomForestClassifier

# Library for Support Vector Machine 
from sklearn.svm import SVC

#Suppress lots of memory warning
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 
from os.path import expanduser as eu
from os import path

scanpy.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)

from itertools import chain

##################################
# Define command line parameters #
##################################

parser = argparse.ArgumentParser()

# outdir
parser.add_argument('-o', '--outdir', nargs='+',  type=str, help='Path where the supervised model predictions on the target dataset are saved.', required = True)

# from
parser.add_argument('-f', '--adata-from', type=str, nargs='+', help='Path to dataset we want to transfer labels from. This must be an anndata object with rows as cells and genes as columns.',required=True)

# to
parser.add_argument('-t','--adata-to', type=str, nargs='+', help='Path to dataset we want to transfer labels to. This must be an anndata object with rows as cells and genes as columns.',required=True)

# model
parser.add_argument('-m', '--model', type=str, nargs='+', help='Supervised learning algorithm. Choices are: Support Vector Machine (SVM), Logistic Regression (LR), Random Forest (RF).', choices=['SVM', 'LR', 'RF'], required=True)

# is-raw
parser.add_argument('-ir', '--is-raw', type=str, nargs='+', help='Boolean value specifying if the anndatas in --from and --to are raw data. If True, the program will skip the resetting to raw step.',default=True)

# n-hvgs
parser.add_argument('-n', '--n-hvgs', type=int, help='Number of highly variable genes to use as features to train the classifier.',required = True)

# downsample (this is mostly for computational efficiency)
parser.add_argument('-d', '--downsample', type=int, help='Maximum number of instances per class. If a class has more instances than the specified number, it will be downsampled.', default=0)

# weights (adjust class weights for unbalanced dataset)
parser.add_argument('-w', '--weights', type=str, nargs='+', help='Boolean value specifying whether or not to run the supervised learning algorithm with class weights. If True, applies smoothing weights to increase the power of minority classes and reduce power of majority classes.', default=True)

# labels
parser.add_argument('-l', '--labels', type=str, nargs='+', help='Name of the column in the metadata of the anndata specified in --from that contains the labels we wish to transfer to the anndata in --to.', required=True)

args = parser.parse_args()


####################
# Print parameters #
####################

print("Launching label transfer script with arguments:/n {}".format(args))

print('''
#################
# Load datasets #
#################
''')


# Load dataset we want to transfer label from
adata_from = scanpy.read(args.adata_from[0])
print("adata_from dims: {}".format(adata_from.shape))

# Load dataset we want to transfer label to
adata_to = scanpy.read(args.adata_to[0])
print("adata_to dims: {}".format(adata_to.shape))   

# Re-set to raw data if needed 
def reset_raw(adata):
    adata = anndata.AnnData(X = adata.raw.X, var = adata.raw.var, obs = adata.obs)
    return adata

if args.is_raw[0] != 'True':
    adata_from = reset_raw(adata_from)
    adata_to = reset_raw(adata_to)
    print('Datasets have been reset to raw counts.')
else: 
    print('Datasets have not been reset to raw counts as the user has not specified the need for it.')


##############
# Downsample #
##############

# Downsample for computational efficiency (if needed)
def downsample(adata, labels): 
    
    myindex = adata.obs[labels].value_counts().index 
    myvalues = adata.obs[labels].value_counts().values
    clusters = pandas.Series(myvalues, index = myindex)
    
    # Find clusters with > n cells (n is user-defined and specified by the argument --downsample)
    n = eu(args.downsample)
    cl2downsample = clusters.index[ clusters.values > n ]

    # save all barcode ids from small clusters
    holder = []
    holder.append( adata.obs_names[[ i not in cl2downsample for i in adata.obs[labels] ]] ) 

    # randomly sample n cells in the cl2downsample
    for cl in cl2downsample:
        print(cl)
        cl_sample = adata[[ i == cl for i in adata.obs[labels]]].obs_names
        # n = int(round(len(cl_sample)/2, 0))
        cl_downsample = random.sample(set(cl_sample), n )
        holder.append(cl_downsample)
    
    # samples to include
    samples = list(chain(*holder))

    # Filter adata_count
    adata = adata[[ i in samples for i in adata.obs_names ]]
    return adata

if args.downsample > 0:
    adata_from = downsample(adata_from, args.downsample)
    print("Dataset has been downsampled")
else: 
    print('Dataset has not been downsampled as the user has not specified the need for it.')
        

print('''
###################
# Intersect genes #
###################
''') 

adata_from_genes = adata_from.var_names.to_list()
adata_to_genes = adata_to.var_names.to_list()

from functools import reduce
inters = reduce(numpy.intersect1d, (adata_from_genes, adata_to_genes))
print('There are {} genes that are shared between the two datasets (before removing cell cycle-associated genes)'.format(len(inters)))

print('''
###########################
# Remove cell cycle genes #
###########################
''') 

cell_cycle_genes = [x.strip() for x in open(file='/nfs/users/nfs_v/vl6/regev_lab_cell_cycle_genes.txt')]
cell_cycle_genes = [x for x in cell_cycle_genes if x in list(inters)]
inters = [x for x in list(inters) if x not in cell_cycle_genes]
print('There are {} genes that are shared between the two datasets (after removing cell cycle-associated genes)'.format(len(inters)))

adata_from = adata_from[:, list(inters)]

print('''
######################
# Data preprocessing #
######################
''') 

def preprocessing_training(adata, hvgs):
    # Per cell normalization
    scanpy.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    # Log transformation 
    scanpy.pp.log1p(adata)
    # Filter HVGs --> Select top N highly variable genes that will serve as features to the machine learning models  
    scanpy.pp.highly_variable_genes(adata, n_top_genes = hvgs)
    highly_variable_genes = adata.var["highly_variable"]
    adata = adata[:, highly_variable_genes]
    # Scale
    scanpy.pp.scale(adata, max_value=10)
    return adata

def preprocessing_testing(adata, hvgs): 
    # Per cell normalization
    scanpy.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    # Log transformation 
    scanpy.pp.log1p(adata)
    # Use the HVGs that were selected for the training set 
    adata = adata[:, hvgs]
    # Scale
    scanpy.pp.scale(adata, max_value=10)
    return adata

X_df = pandas.DataFrame(adata_from.X.toarray(), index = adata_from.obs_names, columns = adata_from.var_names)  # Fetching the count matrix to use as input to the model 
print(type(X_df), X_df.shape)

# Labels we want to predict 
y = list(adata_from.obs[args.labels[0]].astype('str'))
print('List of unique labels we wish to transfer: {}'.format(numpy.unique(y)))

# Split the training dataset into train and test sets 
X_train, X_test, y_train, y_test = train_test_split(
        X_df,
        y,
        test_size=0.25, # This can be changed, though it makes sense to use 25-30% of the data for test
        random_state=1234,
    )

X_train_indices = X_train.index.to_list()
X_test_indices = X_test.index.to_list()
X_train = preprocessing_training(adata_from[X_train_indices, :], args.n_hvgs)

hvgs_X_train = X_train.var_names.to_list()
X_test = preprocessing_testing(adata_from[X_test_indices, :], hvgs_X_train)

X_train = X_train.X
X_test = X_test.X

print('''
####################################
# Train supervised learning models #
####################################
''') 

# Compute smoothing class weights 
def class_weight(labels_dict,mu=0.15):
    total = numpy.sum(list(labels_dict.values()))
    keys = labels_dict.keys()
    weight = {}
    print(total)

    for i in keys:
        score = numpy.log(mu*total/float(labels_dict[i]))
        weight[i] = score if score > 1 else 1
    return weight

labels_dict = {}
for i in y_train:
    labels_dict[i] = labels_dict.get(i, 0) + 1

print(labels_dict)
weights = class_weight(labels_dict)

if args.model[0] == 'SVM':
    print('Training a Support Vector Machine multiclass classifier')
    # Instantiate an RBF Support Vector Machine
    if args.weights == True:
        print("Applying smoothing class weights")
        svm = SVC(kernel = "rbf", probability = True, class_weight = weights)
    else:
        print("NOT applying smoothing class weights")
        svm = SVC(kernel = "rbf", probability = True, class_weight = None)

    # Instantiate a PCA 
    pca = PCA()

    # Create pipeline object
    pipe = Pipeline(steps=[('pca', pca), ('SVC', svm)])

    print('Starting hyperparameter tuning with exhaustive grid search')

    # Choose a grid of hyperparameters values (these are arbitrary but reasonable as I took reference values from the documentation)
    params_svm = {'SVC__C':[0.1, 1, 10, 100], 'SVC__gamma':[0.001, 0.01, 0.1], 'pca__n_components': [0.7, 0.8, 0.9]}

    # Use grid search cross validation to span the hyperparameter space and choose the best 
    grid = GridSearchCV(pipe, param_grid = params_svm, cv=5, verbose =1, n_jobs = -1)

    # Fit the model to the training set of the training data
    grid.fit(X_train, y_train)

    # Report the best hyperparameters and the corresponding score
    print("Best Support Vector Machine CV params", grid.best_params_)
    print("Best Support Vectore Machine CV accuracy", grid.best_score_)
    
elif args.model[0] == 'LR':
    print('Training a Logistic Regression multiclass classifier')
    # Instantiate a Logistic Regression Classifier and specify L2 regularization
    if args.weights == True:
        print("Applying smoothing class weights")
        lr = LogisticRegression(penalty='l2', multi_class="multinomial", max_iter = 2000, class_weight = weights)
    else: 
        lr = LogisticRegression(penalty='l2', multi_class="multinomial", max_iter = 2000, class_weight = None)

    # Instantiate a PCA object
    pca = PCA()

    # Create pipeline object
    pipe = Pipeline(steps=[('pca', pca), ('LogReg', lr)])

    print('Hyperparameter tuning with randomized grid search')

    # Choose a grid of hyperparameters values (these are arbitrary but reasonable as I took reference values from the documentation)
    params_lr = {'LogReg__C' : [0.001, 0.01, 0.1, 1, 10, 100], 'LogReg__solver' : ["lbfgs", 'newton-cg', 'sag'], 
                               'pca__n_components' : [0.7, 0.8, 0.9]}

    # Use grid search cross validation to span the hyperparameter space and choose the best 
    grid = GridSearchCV(estimator = pipe, param_grid =  params_lr, cv = 5, n_jobs = -1)

    # Fit the model to the training set of the training data
    grid.fit(X_train, y_train)

    # Report the best parameters
    print("Best Logistic Regression CV params", grid.best_params_)

    # Report the best hyperparameters and the corresponding score
    print("Softmax train accuracy:", grid.score(X_train, y_train))
    print("Softmax test accuracy:", grid.score(X_test, y_test))

elif args.model[0] == 'RF':
    print('Training a Random Forest multiclass classifier')
    # Instantiate a Random Forest Classifier 
    if args.weights == True:
        print("Applying smoothing class weights")
        rf = RandomForestClassifier(class_weight = weights)
    else: 
        rf = RandomForestClassifier(class_weight = None)
    # Instantiate a PCA object
    pca = PCA()
    # Create pipeline object
    pipe = Pipeline(steps=[('pca', pca), ('RF', rf)])

    print('Hyperparameter tuning with randomized grid search')

    # Choose a grid of hyperparameters values (these are arbitrary but reasonable as I took reference values from the documentation)
    params_rf = {"RF__n_estimators": [50, 100, 200, 300], 'RF__min_samples_leaf': [1, 5], 'RF__min_samples_split': [2, 5, 10], 
                 'pca__n_components' : [0.7, 0.8,0.9]}

    # Use grid search cross validation to span the hyperparameter space and choose the best 
    grid = GridSearchCV(estimator = pipe, param_grid = params_rf, cv = 5, n_jobs = -1)

    # Fit the model to the training set of the training data
    grid.fit(X_train, y_train)
    
    # Report the best hyperparameters and the corresponding score
    print("Best Random Forest CV params", grid.best_params_)
    print("Best Random Forest CV accuracy", grid.best_score_)
else: 
    print('The model you specified in either invalid or not supported by this script!')
    
            
        
print('''
#######################################
# Evaluate supervised learning models #
#######################################
''') 

predicted_labels = grid.best_estimator_.predict(X_test) 
report = classification_report(y_test, predicted_labels)
print("Accuracy:", accuracy_score(y_test, predicted_labels))
print(report)

# Confusion matrix 
cnf_matrix = confusion_matrix(y_test, predicted_labels)

class_names=[0,1] # name  of classes
fig, ax = matplotlib.pyplot.subplots()
tick_marks = numpy.arange(len(class_names))
matplotlib.pyplot.xticks(tick_marks, class_names)
matplotlib.pyplot.yticks(tick_marks, class_names)
# create heatmap
seaborn.heatmap(pandas.DataFrame(cnf_matrix), annot=True, cmap="YlGnBu" ,fmt='g')
ax.xaxis.set_label_position("top")
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.title('Confusion matrix', y=1.1)
matplotlib.pyplot.ylabel('Actual label')
matplotlib.pyplot.xlabel('Predicted label')
matplotlib.pyplot.savefig(args.outdir[0] + 'confusion_matrix.png', bbox_inches='tight')


print('''
###################
# Transfer labels #
###################
''')

grid.best_estimator_.feature_names = hvgs_X_train

def process_and_subset_data(adata, genes):
    # save the log transformed counts as raw 
    adata.raw = adata.copy()
    # Per cell normalization
    scanpy.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    # Log transformation 
    scanpy.pp.log1p(adata)
    # Subset data
    adata = adata[:, list(genes)]
    # Scale
    scanpy.pp.scale(adata, max_value=10)
    return adata

def make_single_predictions(adata, classifier): 
    #if scipy.sparse.issparse(adata.X):
        #adata.X = adata.X.toarray()
    adata_X = numpy.array(adata.X)
    print(type(adata_X), adata_X.shape)
    adata_preds = classifier.predict(adata_X)
    adata.obs['classifier'] = adata_preds
    print(adata.obs.classifier.value_counts(dropna = False))
    
def make_correspondence(classifier):
    corr = {}
    for i in range(0,len(classifier.classes_)):
            corr[i] = classifier.classes_[i]
    return corr

def make_probability_predictions(adata, classifier):
    adata_X = numpy.array(adata.X)
    print(type(adata_X), adata_X.shape)
    proba_preds = classifier.predict_proba(adata_X)
    df_probs = pandas.DataFrame(numpy.column_stack(list(zip(*proba_preds))))
    corr = make_correspondence(classifier)
    for index in df_probs.columns.values:
        celltype = corr[index]
        adata.obs['prob_'+celltype] = df_probs[index].to_list()

adata_to = process_and_subset_data(adata_to, grid.best_estimator_.feature_names)
make_single_predictions(adata_to, grid.best_estimator_)
make_probability_predictions(adata_to, grid.best_estimator_)

to_save = ['classifier']
for i in adata_to.obs.columns:
    if i.startswith('prob_'):
        to_save.append(i)

print('Saving the following obs: {}'.format(to_save))

print('''
###########################
# Save transferred labels #
###########################
''')

adata_to.obs[to_save].to_csv(args.outdir[0] + 'transferred_labels.csv')
