
"""
Given a gene list (rnk file) and list of gene sets (gmt file), perform GSPA to identify significant gene sets.

Example command-line usage:
$ python gspa.py --rnk_file rnk_files/GSE4183.rnk --gmt_file gene_sets/kegg.gmt --output_folder output_folder

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

"""

from functions import *
# import pickle
# import pandas as pd
# import numpy as np
# import argparse

if __name__ == '__main__':

    #%% Load command-line arguments
    parser = argparse.ArgumentParser(description='Run GSPA')

    parser.add_argument('--rnk_file', type=str, required=True, help='Path to rnk file')
    parser.add_argument('--gmt_file', type=str, required=True, help='Path to gmt file')
    parser.add_argument('--output_folder', type=str, required=False, default='output_folder', help = 'Directory to store output')
    parser.add_argument('--use_permutation', type=bool, required = False, default=True, 
                        help='For generating null distributions, whether to use permutation of the seed set ("True", default), or permutation of the expanded set ("False").')
    parser.add_argument('--r', type=float, required=False, default=0.3, help='Radius for set expansion (default = 0.3)')
    parser.add_argument('--min_set_size', type=int, required=False, default=2, help='Minimum gene set size to consider (default = 2)')
    parser.add_argument('--max_set_size', type=int, required=False, default=300, help='Maximum gene set size to consider (default = 300)')
    parser.add_argument('--n_perm', type=int, required=False, default=100, help='Number of permutations for generating null distributions (default = 100)')
    parser.add_argument('--weighted_score_type', type=int, required=False, default=1, help='Exponential degree of weighting of each element for weighted K-S test (default = 1)')
    parser.add_argument('--verbose', type=bool, required=False, default=True, help='Whether to print progress (default = True)')
    parser.add_argument('--drop_na', type=bool, required=False, default=True, help='If True (default), drop gene sets for which a p-value could not be computed.')

    args = parser.parse_args()

    rnk_file = args.rnk_file
    gmt_file = args.gmt_file
    output_folder = args.output_folder
    radius = args.r
    min_set_size = args.min_set_size
    max_set_size = args.max_set_size
    nperm = args.n_perm
    weighted_score_type = args.weighted_score_type
    use_permutation = args.use_permutation
    verbose = args.verbose
    drop_na = args.drop_na

    #%% Load embeddings
    node_ids = pickle.load(open('embeddings/humanppi_node_ids.p', 'rb'))
    node_embeddings = pickle.load(open('embeddings/humanppi_node_embs.p', 'rb'))

    #%% Load gene sets
    gene_set_dict = gmt2dict(gmt_file)

    #%% Make output directory
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    #%% Load gene list
    gene_list_name = rnk_file.split('/')[-1].split('.')[-2]
    gene_list_input = list(pd.read_csv(rnk_file,sep='\t',header=None)[0])
    correl_vector_input = list(pd.read_csv(rnk_file,sep='\t',header=None)[1])

    print(f'{gene_list_name} loaded.')

    #%%
    results = []

    print('Running GSPA...')

    # Check which user-supplied genes have embeddings
    gene_list_innetwork_indices = np.arange(len(gene_list_input))[np.in1d(gene_list_input, node_ids, assume_unique=True)]
    gene_list_outnetwork_indices = np.arange(len(gene_list_input))[~np.in1d(gene_list_input, node_ids, assume_unique=True)]
    gene_list_innetwork = np.take(gene_list_input, gene_list_innetwork_indices)

    # Precompute all gene-gene embedding distances
    XA_full = np.array([node_embeddings[node_ids.index(x)] for x in gene_list_innetwork])
    d_full = distance.cdist(XA_full, XA_full, 'cosine')

    XA_big = np.array([node_embeddings[node_ids.index(x)] for x in node_ids])
    XB_big = np.array([node_embeddings[node_ids.index(x)] for x in gene_list_innetwork])
    d_big = distance.cdist(XA_big, XB_big, 'cosine')

    # For each gene set, run GSPA algorithm to generate ES, NES, P-value, and FDR
    nsets = len(gene_set_dict.keys())
    for i,gene_set_name in enumerate(gene_set_dict.keys()):
        try:
            if gene_set_name not in [x[0] for x in results]:
                if i < 10000: # maximum number of sets to test
                    gene_set_input = gene_set_dict[gene_set_name]
                    gene_list_input = gene_list_input
                    gs = [x for x in gene_set_input if x in node_ids]
                    if (len(gs)>min_set_size) & (len(gs)<max_set_size):
                        es_0, esnull_0, hit_ind_0, RES_0, gsl_0 = enrichment_score_gspa(gene_list_input, correl_vector_input, gene_set_input, 
                                                        node_ids, node_embeddings, radius=radius,
                                                        weighted_score_type=weighted_score_type, nperm=nperm, 
                                                        rs=None, single=False, scale=False, plot=False,
                                                        use_permutation=use_permutation,
                                                        d_big = d_big, 
                                                        d_full = d_full,
                                                        gene_list_innetwork_indices=gene_list_innetwork_indices, 
                                                        gene_list_outnetwork_indices=gene_list_outnetwork_indices,
                                                        gene_list_innetwork = gene_list_innetwork)
                        results.append((
                                gene_set_name,len(gene_set_input),
                                es_0,np.array(esnull_0),np.array(hit_ind_0),np.array(RES_0),
                                ))
                        if verbose:
                            print(i+1, 'of', nsets, gene_set_name)
                    else:
                        # print(gene_set_name + ' is wrong size.')
                        continue
            else:
                break
        except Exception as e:
            print(f'Skipped gene set {gene_set_name} with error {e}')
            pass


    #%% Convert results to dataframe
    name = gene_list_name + '_results'
    es_s = [x[2] for x in results]
    esnull_s = [x[3].tolist() for x in results]
    gs = gspa_significance(es_s, esnull_s)
    df_results = pd.DataFrame([x for x in gs],columns=['ES','NES','P-value','FDR'])
    df_results['gene_set_name'] = pd.DataFrame([x[0] for x in results], index=df_results.index)
    # df_results['size0'] = pd.DataFrame([x[1] for x in results], index=df_results.index)
    # df_results['size1'] = pd.DataFrame([x[10] for x in results], index=df_results.index)

    df_results.rename(columns={'gene_set_name': 'Gene Set'}, inplace=True)

    df_results.set_index('Gene Set',inplace=True,drop=True)

    # es_s = [x[6] for x in results]
    # esnull_s = [x[7].tolist() for x in results]
    # gs = gspa_significance(es_s, esnull_s)
    # df_results[['es_1','nes_1','pval_1','fdr_1']] = pd.DataFrame([x for x in gs], index=df_results.index)
    # # df_results['size'] = pd.DataFrame([x[1] for x in results], index=df_results.index)

    if drop_na:
        df_results.replace(np.inf, np.nan, inplace=True)
        df_results.dropna(inplace=True)

    df_results.sort_values(by='FDR',ascending=True,inplace=True)

    print(f'Job complete. Results saved to {output_folder}.')

    df_results.to_csv(f'{output_folder}/'+name+'.csv')

