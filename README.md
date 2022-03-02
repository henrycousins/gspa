# Gene Set Proximity Analysis (GSPA)

Gene set proximity analysis (GSPA) is a method for identifying critical gene sets in functional genetic datasets using low-dimensional gene embeddings. For additional documentation, please refer to our preprint, available at https://arxiv.org/abs/2202.00143. 

![Alt text](overview.png?raw=true "GSPA Overview")

# Dependencies
GSPA is built in Python 3 and requires Matplotlib, NumPy, pandas, and SciPy, in addition to the standard Python library. It was tested using Python 3.7, Matplotlib 3.3.4, NumPy 1.20.2, pandas 1.2.4, and SciPy 1.6.2.

# Installation

1. Clone the repository

```
$ git clone https://github.com/henrycousins/gspa
$ cd gspa
```

2. Ensure that NumPy, pandas, SciPy, and Matplotlib are installed. We provide a conda environment containing all necessary dependencies, which can be activated by running

```
$ conda env create -f environment.yml
$ conda activate GSPA
```
# Simple usage

GSPA can be run through the command line as follows:

```
python gspa.py --rnk_file RNK_FILE --gmt_file GMT_FILE --output_folder OUTPUT_FOLDER
```
where ```RNK_FILE``` is the path to the .rnk file containing the gene list, ```GMT_FILE``` is the path to the .gmt file containing the gene sets, and ```OUTPUT_FOLDER``` is name of the directory where results should be stored. For instance, the following code performs GSPA on the GSE4183 gene list using KEGG pathway gene sets, saving the results to ```/outputs```:

```
python gspa.py --rnk_file rnk_files/GSE4183.rnk --gmt_file gene_sets/kegg.gmt --output_folder outputs
```

# Advanced usage

Many additional parameters can be defined by the user.

```
usage: gspa.py [-h] --rnk_file RNK_FILE --gmt_file GMT_FILE
               [--output_folder OUTPUT_FOLDER]
               [--use_permutation USE_PERMUTATION] [--r R]
               [--min_set_size MIN_SET_SIZE] [--max_set_size MAX_SET_SIZE]
               [--n_perm N_PERM] [--weighted_score_type WEIGHTED_SCORE_TYPE]
               [--verbose VERBOSE] [--drop_na DROP_NA]

optional arguments:
  -h, --help            show this help message and exit
  --rnk_file RNK_FILE   Path to rnk file
  --gmt_file GMT_FILE   Path to gmt file
  --output_folder OUTPUT_FOLDER
                        Directory to store output
  --use_permutation USE_PERMUTATION
                        For generating null distributions, whether to use
                        permutation of the seed set ("True", default), or
                        permutation of the expanded set ("False").
  --r R                 Radius for set expansion (default = 0.3)
  --min_set_size MIN_SET_SIZE
                        Minimum gene set size to consider (default = 2)
  --max_set_size MAX_SET_SIZE
                        Maximum gene set size to consider (default = 300)
  --n_perm N_PERM       Number of permutations for generating null
                        distributions (default = 100)
  --weighted_score_type WEIGHTED_SCORE_TYPE
                        Exponential degree of weighting of each element for
                        weighted K-S test (default = 1)
  --verbose VERBOSE     Whether to print progress (default = True)
  --drop_na DROP_NA     If True (default), drop gene sets for which a p-value
                        could not be computed.
```

