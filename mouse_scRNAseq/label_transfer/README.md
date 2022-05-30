### 1: Converting single-cell dataset formats 
**Seurat2Anndata.R**

Parameters: 
 * ``--seurat_object`` : Path to Seurat object you wish to convert to Anndata (string)
 * ``--object_name`` : Name of the saved object (string)
 * ``--outdir`` : Path to where you want to store the resulting Anndata object (string) 

### 2: Converting gene IDs across species (mouse to human only for now)
**CrossSpecies_Biomart.R**

Parameters: 
 * ``--genes``: Path to txt file with mouse ENSEMBL gene IDs you wish to convert to human (string)
 * ``--outdir``: Path to where you want to store the resulting converted human ENSEMBL IDs and gene names (string)

### 3: Transfering labels between single-cell datasets with supervised machine learning algorithms 
**LabelTransfer_SupervisedLearning.py** 

Parameters: 
 * ``--adata-from`` (-f) : Path to dataset we want to transfer labels from. This must be an anndata object with rows as cells and genes as columns. (string)
 * ``--adata-to`` (-t) : Path to dataset we want to transfer labels to. This must be an anndata object with rows as cells and genes as columns. (string)
 * ``--model`` (-m) : Supervised learning algorithm. Choices are: SVM (Support Vector Machine), LR (Logistic Regression), RF (Random Forest). (string) 
 * ``--is-raw`` (-ir) : Boolean value specifying if the anndatas in --from and --to are raw data. If True, the program will skip the resetting to raw step. (string)
 * ``--n-hvgs`` (-n) : Number of highly variable genes to use as features to train the classifier. (integer) 
 * ``--downsample`` (-d) : Maximum number of instances per class. If a class has more instances than the specified number, it will be downsampled. Default is 0, meaning no downsampling will be performed. (int) 
 * ``--weights`` (-w) : Boolean value specifying whether or not to run the supervised learning algorithm with class weights. If True, applies smoothing weights to increase the power of minority classes and reduce power of majority classes. (string)
 * ``--labels`` (-l) : Name of the column in the metadata of the anndata specified in --adata-from that contains the labels we wish to transfer to the anndata in --adata-to. (string)
 * ``--outdir`` (-o) : Path where the supervised model predictions on the target dataset are saved. (string)

The script requires the python libraries specified in **requirements.txt**. If you would like to install them use:

``pip install -r requirements.txt``
