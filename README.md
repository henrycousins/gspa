# Gene Set Proximity Analysis (GSPA)

Gene set proximity analysis (GSPA) is a method for identifying critical gene sets in functional genetic datasets using low-dimensional gene embeddings. For additional documentation, please refer to our preprint, available at https://arxiv.org/abs/2202.00143. 

![Alt text](overview.png?raw=true "GSPA Overview")

## Dependencies
GSPA is built in Python 3 and requires Matplotlib, NumPy, pandas, and SciPy, in addition to the standard Python library. It was tested using Python 3.7, Matplotlib 3.3.4, NumPy 1.20.2, pandas 1.2.4, and SciPy 1.6.2. We provide an optional conda environment containing all necessary dependencies, which can be activated by running the following:
```
$ conda env create -f environment.yml
$ conda activate GSPA
```

## Installation

Once dependencies are installed, clone the repository as follows:

```
$ git clone https://github.com/henrycousins/gspa
$ cd gspa
```

## Simple usage

GSPA can be run through the command line as follows:

```
python gspa.py --rnk_file RNK_FILE --gmt_file GMT_FILE --output_folder OUTPUT_FOLDER
```
where ```RNK_FILE``` is the path to the .rnk file containing the gene list, ```GMT_FILE``` is the path to the .gmt file containing the gene sets, and ```OUTPUT_FOLDER``` is name of the directory where results should be stored. (Please note that ```OUTPUT_FOLDER``` is created automatically.) For instance, the following code performs GSPA on the GSE4183 gene list using KEGG pathway gene sets, saving the results to ```/outputs```:

```
python gspa.py --rnk_file rnk_files/GSE4183.rnk --gmt_file gene_sets/kegg.gmt --output_folder outputs
```

## Advanced usage

Additional parameters can be defined by the user.

```
usage: gspa.py [-h] --rnk_file RNK_FILE --gmt_file GMT_FILE
               [--output_folder OUTPUT_FOLDER]
               [--use_permutation USE_PERMUTATION] [--r R]
               [--min_set_size MIN_SET_SIZE] [--max_set_size MAX_SET_SIZE]
               [--max_num_sets MAX_NUM_SETS] [--n_perm N_PERM]
               [--weighted_score_type WEIGHTED_SCORE_TYPE] [--verbose VERBOSE]
               [--drop_na DROP_NA] [--results_file RESULTS_FILE]

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
  --r R                 Radius for set expansion (default = 0.1)
  --min_set_size MIN_SET_SIZE
                        Minimum gene set size to consider (default = 2)
  --max_set_size MAX_SET_SIZE
                        Maximum gene set size to consider (default = 300)
  --max_num_sets MAX_NUM_SETS
                        Maximum number of gene sets to evaluate (default =
                        600)
  --n_perm N_PERM       Number of permutations for generating null
                        distributions (default = 100)
  --weighted_score_type WEIGHTED_SCORE_TYPE
                        Exponential degree of weighting of each element for
                        weighted K-S test (default = 1)
  --verbose VERBOSE     Whether to print progress (default = True)
  --drop_na DROP_NA     If True, drop gene sets for which an exact NES could
                        not be computed (default = False).
  --results_file RESULTS_FILE
                        Custom name for results file (default =
                        [list_name]_results.csv).
```

## Detailed file descriptions

File | Type | Extension | Columns | Comments
--- | --- | --- | --- | ---
```RNK_FILE``` | Input | .rnk | [gene ID, gene score] | [Description](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#RNK:_Ranked_list_file_format_.28.2A.rnk.29)
```GMT_FILE``` | Input | .gmt | [gene set name, description (optional), gene 1, gene 2, ... , gene n] | [Description](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29)
```RESULTS_FILE``` | Output | .csv | [gene set name, ES, NES, P-value, FDR] | Created automatically

```
.
├── backup
└── gspa
    ├── LICENSE           # LICENSE
    ├── README.md         # README
    ├── embeddings        # contains gene embeddings and IDs as pickled arrays
    ├── environment.yml   # optional dependency environment
    ├── functions.py      # helper functions
    ├── gene_sets         # contains example gene sets
    ├── gspa.py           # main
    ├── outputs           # example directory for results
    ├── overview.png      # overview graphic
    ├── rnk_files         # contains example rnk files
    └── data              # supplementary data
        ├── ppi           # related to ppi network and embeddings
        ├── sarscov2      # related to sars-cov-2 screens
        └── set_matching  # related to comparing similar gene sets
```

## License

This software is available under an MIT license.

## Sustainability

We anticipate releasing core updates on a 6-month basis for a minimum of 2 years and will respond to any issues regarding the codebase, which is open-source, as they arise. Supporting funding for GSPA is provided by R01 GM102365 from National Institutes of Health.