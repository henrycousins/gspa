
#%%
import os
# import re
import numpy as np
import pandas as pd
# import glob as glob
# import scipy.stats
import pickle
import random as random
# from sklearn import metrics as skm
# import seaborn as sns
from scipy.spatial import distance
# from sklearn import preprocessing, model_selection
# from IPython.display import display, HTML
import matplotlib.pyplot as plt
import sys, logging
# import numpy as np
#from functools import reduce
#from multiprocessing import Pool
from math import ceil
import argparse
# import matplotlib.pyplot as plt
# from scipy.spatial import distance


#%%


gene_list_input = None
gene_list_innetwork = None
d_full = None
gene_list_innetwork_indices = None
d_big = None
# gene_list_innetwork_indices = None
gene_list_outnetwork_indices = None


def gmt2dict(gmt):
    """
    Load gmt file as dictionary
    """
    with open(gmt) as genesets:
        genesets_dict = { line.strip().split("\t")[0]: line.strip().split("\t")[2:]
                                    for line in genesets.readlines()}
    return genesets_dict

# gene_list_input = None
# gene_list_innetwork = None
# d_full = None
# gene_list_innetwork_indices = None

def compute_permutation(randomstate, radius, tag_indicator_b, gene_list_outnetwork_indices,
                        gene_list_innetwork_indices = gene_list_innetwork_indices,
                        gene_list_input = gene_list_input, 
                        gene_list_innetwork = gene_list_innetwork, 
                        d_full = d_full):
    """
    Generate null distribution by permuting seeds and expanding to proximal genes
    """

    rrs = np.random.RandomState(randomstate)
    tag_indicator = np.zeros(len(gene_list_input))

    if radius == 0:
        iii = rrs.choice(tag_indicator.shape[0],size=sum(tag_indicator_b))
        tag_indicator[iii] = 1

    else:
        infrac = gene_list_innetwork_indices.shape[0] / len(gene_list_input)
        outfrac = gene_list_outnetwork_indices.shape[0] / len(gene_list_input)

        randomgeneset_in = rrs.choice(gene_list_innetwork.shape[0],size=int(sum(tag_indicator_b)*infrac))

        d_in = d_full[randomgeneset_in,:]

        tag_indicator_c = np.amin(d_in,axis=0)

        tag_indicator_c = np.array([1 if d <= radius else 0 for d in tag_indicator_c])
        tag_indicator[gene_list_innetwork_indices] = tag_indicator_c

        randomgeneset_out = rrs.choice(gene_list_outnetwork_indices,size=int(sum(tag_indicator_b)*outfrac))
        tag_indicator[randomgeneset_out] = 1

    return tag_indicator

# d_big = None
# # gene_list_innetwork_indices = None
# gene_list_outnetwork_indices = None

def enrichment_score_gspa(gene_list_input, correl_vector_input, gene_set_input, node_ids, node_embeddings, radius=0, 
            weighted_score_type=1, nperm=1000, rs=None, single=False, scale=False, plot=False, pval=False,
            use_permutation=True, 
            d_big = d_big, 
            d_full = d_full,
            gene_list_innetwork_indices=gene_list_innetwork_indices, 
            gene_list_outnetwork_indices=gene_list_outnetwork_indices,
            gene_list_innetwork = gene_list_innetwork):
    
    """
    Compute ES using the GSPA algorithm
    """

    radius = radius

    gene_set_innetwork = [x for x in gene_set_input if x in node_ids]

    gene_set_innetwork_nodeindices = np.arange(len(node_ids))[np.in1d(node_ids, gene_set_innetwork, assume_unique=True)]

    correl_vector = correl_vector_input

    N = len(gene_list_input)

    tag_indicator = np.zeros(len(gene_list_input))

    tag_indicator_b = np.in1d(gene_list_input, gene_set_input, assume_unique=True).astype(int)

    d = d_big[gene_set_innetwork_nodeindices,:]

    tag_indicator_c = np.amin(d,axis=0)

    tag_indicator_c = np.array([1 if d <= radius else 0 for d in tag_indicator_c])

    tag_indicator[gene_list_innetwork_indices] = tag_indicator_c
    tag_indicator = np.add(tag_indicator,tag_indicator_b)
    tag_indicator = np.where(tag_indicator>0,1,0)

    if weighted_score_type == 0 :
        correl_vector = np.repeat(1, N)
    else:
        correl_vector = np.abs(correl_vector)**weighted_score_type

    # get indices of tag_indicator
    hit_ind = np.flatnonzero(tag_indicator).tolist()
    # if used for compute esnull, set esnull equal to permutation number, e.g. 1000
    # else just compute enrichment scores
    # set axis to 1, because we have 2D array
    axis = 1
    tag_indicator = np.tile(tag_indicator, (nperm+1,1))
    correl_vector = np.tile(correl_vector,(nperm+1,1))
    # gene list permutation
    rs = np.random.RandomState(rs)

    if use_permutation:
        for i in range(nperm):
            tag_indicator[i] = compute_permutation(randomstate=i, radius=radius, 
                                                    tag_indicator_b=tag_indicator_b, 
                                                    gene_list_outnetwork_indices=gene_list_outnetwork_indices,
                                                    gene_list_innetwork_indices=gene_list_innetwork_indices,
                                                    gene_list_input = gene_list_input, 
                                                    gene_list_innetwork = gene_list_innetwork, 
                                                    d_full = d_full)
    else:
        for i in range(nperm): rs.shuffle(tag_indicator[i])

    # for i in range(nperm): rs.shuffle(tag_indicator[i])
    # np.apply_along_axis(rs.shuffle, 1, tag_indicator)

    Nhint = tag_indicator.sum(axis=axis, keepdims=True)
    sum_correl_tag = np.sum(correl_vector*tag_indicator, axis=axis, keepdims=True)
    # compute ES score, the code below is identical to gsea enrichment_score method.
    no_tag_indicator = 1 - tag_indicator
    Nmiss =  N - Nhint
    norm_tag =  1.0/sum_correl_tag
    norm_no_tag = 1.0/Nmiss

    RES = np.cumsum(tag_indicator * correl_vector * norm_tag - no_tag_indicator * norm_no_tag, axis=axis)
    #               # indices of overlap  

    if scale: RES = RES / N
    if single:
        es_vec = RES.sum(axis=axis)
    else:
        max_ES, min_ES =  RES.max(axis=axis), RES.min(axis=axis)
        es_vec = np.where(np.abs(max_ES) > np.abs(min_ES), max_ES, min_ES)
    # extract values
    es, esnull, RES = es_vec[-1], es_vec[:-1], RES[-1,:]



    # ###
    # gene_list = list(sorted_data.index)
    # correl_vector = list(sorted_data.t_statistic)
    # gene_set = chr6q21_geneset # supposed to be most enriched set

    # es, esnull, hit_ind, RES = enrichment_score_gspa(gene_list, correl_vector, gene_set, weighted_score_type=1, 
    #                                             nperm=1000, rs=None, single=False, scale=False)


    if plot:

        fig,(ax1,ax2,ax3) = plt.subplots(3,1,sharex=True)
        fig.set_figheight(6)
        fig.set_figwidth(8)
        ax1.vlines(hit_ind,ymin=0,ymax=1,color='k',alpha=0.3)
        ax1.set_title('Gene set')
        ax2.plot([x for x in range(len(correl_vector_input))], correl_vector_input, color='k')
        ax2.set_title('Correlation with phenotype')
        ax3.plot([x for x in range(len(RES))], RES, color='k')
        ax3.set_title('Enrichment score (ES)')
        plt.xlabel('Gene List Rank')
        plt.tight_layout()

        # print(scipy.stats.ttest_1samp(esnull,es))
        print('P-value:\n')
        print(gspa_pval(es,esnull))
        # print(pval_2(es,esnull))

        plt.figure()
        plt.hist(tag_indicator_c, bins=40, alpha=0.5, density = False)
        # plt.hist(tag_indicator_old_distances, bins=10, alpha=0.5, density = False)
        plt.show()

        plt.figure()
        a = plt.hist(esnull,bins=40,density=False)
        
        plt.vlines(es,0,max(a[0]))
        plt.title(es)


    return es, esnull, hit_ind, RES, sum(tag_indicator_c)


def normalize(es, esnull):
    """
    Normalize ES by rescaling positive and negative scores by the mean of ESnull
    """

    nEnrichmentScores =np.zeros(es.shape)
    nEnrichmentNulls=np.zeros(esnull.shape)

    esnull_pos = np.ma.MaskedArray(esnull, mask=(esnull<0)).mean(axis=1)
    esnull_neg = np.ma.MaskedArray(esnull, mask=(esnull>=0)).mean(axis=1)
    esnull_pos = np.array(esnull_pos)
    esnull_neg = np.array(esnull_neg)

    # NES
    nEnrichmentScores  = np.where(es>=0, es/esnull_pos, -es/esnull_neg)
    # NES_NULL
    nEnrichmentNulls = np.where(esnull>=0, esnull/esnull_pos[:,np.newaxis],
                                          -esnull/esnull_neg[:,np.newaxis])
    
    return nEnrichmentScores, nEnrichmentNulls

def gspa_pval(es, esnull):
    """
    Compute empirical p-value
    """

    condlist = [ es < 0, es >=0]
    choicelist = [(esnull < es.reshape(len(es),1)).sum(axis=1)/ (esnull < 0).sum(axis=1),
                  (esnull >= es.reshape(len(es),1)).sum(axis=1)/ (esnull >= 0).sum(axis=1)]
    pvals = np.select(condlist, choicelist)

    return pvals


def gspa_fdr(nEnrichmentScores, nEnrichmentNulls):
    """
    Compute empirical FDR as the fraction of NESnull greater in magnitude than the observed NES
    """

    nvals = np.sort(nEnrichmentNulls.flatten())
    nnes = np.sort(nEnrichmentScores)
    fdrs = []
    for i in range(len(nEnrichmentScores)):
        nes = nEnrichmentScores[i]
        if nes >= 0:
            allPos = int(len(nvals) - np.searchsorted(nvals, 0, side="left"))
            allHigherAndPos = int(len(nvals) - np.searchsorted(nvals, nes, side="left"))
            nesPos = len(nnes) - int(np.searchsorted(nnes, 0, side="left"))
            nesHigherAndPos = len(nnes) - int(np.searchsorted(nnes, nes, side="left"))
        else:
            allPos = int(np.searchsorted(nvals, 0, side="left"))
            allHigherAndPos = int(np.searchsorted(nvals, nes, side="right"))
            nesPos = int(np.searchsorted(nnes, 0, side="left"))
            nesHigherAndPos = int(np.searchsorted(nnes, nes, side="right"))
        try:
            pi_norm = allHigherAndPos / float(allPos)
            pi_obs = nesHigherAndPos / float(nesPos)
            fdr = pi_norm / pi_obs
            fdrs.append(fdr if fdr < 1 else 1.0)
        except:
            fdrs.append(1000000000.0)

    logging.debug("Statistical testing finished.............................")

    return fdrs


def gspa_significance(enrichment_scores, enrichment_nulls):
    """
    Given ES and ESnull, compute NES, P-values, and FDRs
    """

    np.seterr(divide='ignore', invalid='ignore')

    es = np.array(enrichment_scores)
    esnull = np.array(enrichment_nulls)
    logging.debug("Start to compute pvals..................................")

    pvals = gspa_pval(es, esnull).tolist()

    logging.debug("Start to compute nes and nesnull........................")

    nEnrichmentScores, nEnrichmentNulls = normalize(es, esnull)

    logging.debug("Start to compute fdrs..................................")

    fdrs = gspa_fdr(nEnrichmentScores, nEnrichmentNulls)

    return zip(enrichment_scores, nEnrichmentScores, pvals, fdrs)


# # d_big = None
# # # gene_list_innetwork_indices = None
# # gene_list_outnetwork_indices = None

# def enrichment_score_gspa(gene_list_input, correl_vector_input, gene_set_input, node_ids, node_embeddings, radius=0, 
#             weighted_score_type=1, nperm=1000, rs=None, single=False, scale=False, plot=False, pval=False,
#             use_permutation=True, 
#             d_big = d_big, 
#             d_full = d_full,
#             gene_list_innetwork_indices=gene_list_innetwork_indices, 
#             gene_list_outnetwork_indices=gene_list_outnetwork_indices,
#             gene_list_innetwork = gene_list_innetwork):
    
#     """
#     Compute ES using the GSPA algorithm
#     """

#     radius = radius

#     gene_set_innetwork = [x for x in gene_set_input if x in node_ids]

#     gene_set_innetwork_nodeindices = np.arange(len(node_ids))[np.in1d(node_ids, gene_set_innetwork, assume_unique=True)]

#     correl_vector = correl_vector_input

#     N = len(gene_list_input)

#     tag_indicator = np.zeros(len(gene_list_input))

#     tag_indicator_b = np.in1d(gene_list_input, gene_set_input, assume_unique=True).astype(int)

#     d = d_big[gene_set_innetwork_nodeindices,:]

#     tag_indicator_c = np.amin(d,axis=0)

#     tag_indicator_c = np.array([1 if d <= radius else 0 for d in tag_indicator_c])

#     tag_indicator[gene_list_innetwork_indices] = tag_indicator_c
#     tag_indicator = np.add(tag_indicator,tag_indicator_b)
#     tag_indicator = np.where(tag_indicator>0,1,0)

#     if weighted_score_type == 0 :
#         correl_vector = np.repeat(1, N)
#     else:
#         correl_vector = np.abs(correl_vector)**weighted_score_type

#     # get indices of tag_indicator
#     hit_ind = np.flatnonzero(tag_indicator).tolist()
#     # if used for compute esnull, set esnull equal to permutation number, e.g. 1000
#     # else just compute enrichment scores
#     # set axis to 1, because we have 2D array
#     axis = 1
#     tag_indicator = np.tile(tag_indicator, (nperm+1,1))
#     correl_vector = np.tile(correl_vector,(nperm+1,1))
#     # gene list permutation
#     rs = np.random.RandomState(rs)

#     if use_permutation:
#         for i in range(nperm):
#             tag_indicator[i] = compute_permutation(randomstate=i, radius=radius, 
#                                                     tag_indicator_b=tag_indicator_b, 
#                                                     gene_list_outnetwork_indices=gene_list_outnetwork_indices,
#                                                     gene_list_innetwork_indices=gene_list_innetwork_indices,
#                                                     gene_list_input = gene_list_input, 
#                                                     gene_list_innetwork = gene_list_innetwork, 
#                                                     d_full = d_full)
#     else:
#         for i in range(nperm): rs.shuffle(tag_indicator[i])

#     # for i in range(nperm): rs.shuffle(tag_indicator[i])
#     # np.apply_along_axis(rs.shuffle, 1, tag_indicator)

#     Nhint = tag_indicator.sum(axis=axis, keepdims=True)
#     sum_correl_tag = np.sum(correl_vector*tag_indicator, axis=axis, keepdims=True)
#     # compute ES score, the code below is identical to gsea enrichment_score method.
#     no_tag_indicator = 1 - tag_indicator
#     Nmiss =  N - Nhint
#     norm_tag =  1.0/sum_correl_tag
#     norm_no_tag = 1.0/Nmiss

#     RES = np.cumsum(tag_indicator * correl_vector * norm_tag - no_tag_indicator * norm_no_tag, axis=axis)
#     #               # indices of overlap  

#     if scale: RES = RES / N
#     if single:
#         es_vec = RES.sum(axis=axis)
#     else:
#         max_ES, min_ES =  RES.max(axis=axis), RES.min(axis=axis)
#         es_vec = np.where(np.abs(max_ES) > np.abs(min_ES), max_ES, min_ES)
#     # extract values
#     es, esnull, RES = es_vec[-1], es_vec[:-1], RES[-1,:]



#     # ###
#     # gene_list = list(sorted_data.index)
#     # correl_vector = list(sorted_data.t_statistic)
#     # gene_set = chr6q21_geneset # supposed to be most enriched set

#     # es, esnull, hit_ind, RES = enrichment_score_gspa(gene_list, correl_vector, gene_set, weighted_score_type=1, 
#     #                                             nperm=1000, rs=None, single=False, scale=False)


#     if plot:

#         fig,(ax1,ax2,ax3) = plt.subplots(3,1,sharex=True)
#         fig.set_figheight(6)
#         fig.set_figwidth(8)
#         ax1.vlines(hit_ind,ymin=0,ymax=1,color='k',alpha=0.3)
#         ax1.set_title('Gene set')
#         ax2.plot([x for x in range(len(correl_vector_input))], correl_vector_input, color='k')
#         ax2.set_title('Correlation with phenotype')
#         ax3.plot([x for x in range(len(RES))], RES, color='k')
#         ax3.set_title('Enrichment score (ES)')
#         plt.xlabel('Gene List Rank')
#         plt.tight_layout()

#         # print(scipy.stats.ttest_1samp(esnull,es))
#         print('P-value:\n')
#         print(gspa_pval(es,esnull))
#         # print(pval_2(es,esnull))

#         plt.figure()
#         plt.hist(tag_indicator_c, bins=40, alpha=0.5, density = False)
#         # plt.hist(tag_indicator_old_distances, bins=10, alpha=0.5, density = False)
#         plt.show()

#         plt.figure()
#         a = plt.hist(esnull,bins=40,density=False)
        
#         plt.vlines(es,0,max(a[0]))
#         plt.title(es)


#     return es, esnull, hit_ind, RES, sum(tag_indicator_c)




# #%%
# #%%

# node_ids = pickle.load(open('embeddings/humanppi_node_ids.p', 'rb'))
# node_embeddings = pickle.load(open('embeddings/humanppi_node_embs.p', 'rb'))


# #%%

# # gene_set_dict = gmt2dict('gene_sets/kegg_map_2.txt')
# # gene_set_dict = gmt2dict('gene_sets/c2_symbols_gmt.txt')
# gene_set_dict = gmt2dict('gene_sets/D1_gmt.txt')



# # list_of_gene_lists = ['GSE3467','GSE3585','GSE3678']
# # list_of_gene_lists = ['GSE3467']

# # list_of_gene_lists = ['zhu_1','zhu_2','wang_2','wang_3']

# # list_of_gene_lists = ['GSE4107','GSE4183','GSE5281_EC','GSE5281_HIP','GSE5281_VCX']
# # list_of_gene_lists = ['GSE6956AA','GSE6956C','GSE7305','GSE8671','GSE8762','GSE9348','GSE9476','GSE11906','GSE14762']

# # list_of_gene_lists = ['zhu_1', 'zhu_2', 'wang_1','wang_2','wang_3']

# list_of_gene_lists = ['wang_lfc']

# use_absolute_enrichment = False

# # G = pickle.load(open('/embeddings/string_900_nx.p'))

# #%%


# for gene_list_name in list_of_gene_lists:

#     gene_list_input = list(pd.read_csv('rnk_files/'+gene_list_name+'.rnk',sep='\t',header=None)[0])
#     correl_vector_input = list(pd.read_csv('rnk_files/'+gene_list_name+'.rnk',sep='\t',header=None)[1])
    
#     if use_absolute_enrichment==True:

#         rnk = pd.read_csv('rnk_files/'+gene_list_name+'.rnk',sep='\t',header=None)
#         rnk['abs'] = rnk.apply(lambda x: abs(x[1]), axis=1)
#         rnk.sort_values(by='abs',ascending=False,inplace=True)

#         gene_list_input = list(rnk[0])
#         correl_vector_input = list(rnk['abs'])
    
    
    
#     results = []

#     gene_list_innetwork_indices = np.arange(len(gene_list_input))[np.in1d(gene_list_input, node_ids, assume_unique=True)]
#     gene_list_outnetwork_indices = np.arange(len(gene_list_input))[~np.in1d(gene_list_input, node_ids, assume_unique=True)]
#     gene_list_innetwork = np.take(gene_list_input, gene_list_innetwork_indices)

#     XA_full = np.array([node_embeddings[node_ids.index(x)] for x in gene_list_innetwork])
#     d_full = distance.cdist(XA_full, XA_full, 'cosine')


#     XA_big = np.array([node_embeddings[node_ids.index(x)] for x in node_ids])
#     XB_big = np.array([node_embeddings[node_ids.index(x)] for x in gene_list_innetwork])
#     d_big = distance.cdist(XA_big, XB_big, 'cosine')


#     # radius=0.2
#     # ngsea_rerank(gene_list_name)



#     for i,gene_set_name in enumerate(gene_set_dict.keys()):
#         try:
#             if gene_set_name not in [x[0] for x in results]:
#                 if i < 600:
#                     gene_set_input = gene_set_dict[gene_set_name]
#                     gs = [x for x in gene_set_input if x in node_ids]

#                     if (len(gs)>2) & (len(gs)<300):
#                         es_0, esnull_0, hit_ind_0, RES_0, gsl_0 = enrichment_score_gspa(gene_list_input, correl_vector_input, gene_set_input, 
#                                                         node_ids, node_embeddings, radius=0.1,
#                                                         weighted_score_type=1, nperm=100, 
#                                                         rs=None, single=False, scale=False, plot=False,
#                                                         use_permutation=True)
#                         es_1, esnull_1, hit_ind_1, RES_1, gsl_1 = enrichment_score_gspa(gene_list_input, correl_vector_input, gene_set_input, 
#                                                         node_ids, node_embeddings, radius=0.3,
#                                                         weighted_score_type=1, nperm=100, 
#                                                         rs=None, single=False, scale=False, plot=False,
#                                                         use_permutation=True)

#                         results.append(( # changed list to tuple
#                                 gene_set_name,len(gene_set_input), # changed sublists to arrays
#                                 es_0,np.array(esnull_0),np.array(hit_ind_0),np.array(RES_0),
#                                 es_1,np.array(esnull_1),np.array(hit_ind_1),np.array(RES_1),
#                                 gsl_1
#                         ))

#                 # if i % 10 == 0:
#                 #   pickle.dump(results,open(results_name,'wb'))

#                         print(i,gene_list_name, gene_set_name)
#                     else:
#                         print(gene_set_name + ' is wrong size.')
#                         continue
#             else:
#                 break
#         except:
#             print(f'Skipped gene set {gene_set_name}')
#             pass

#     #%%

#     name = gene_list_name + '_results'

#     es_s = [x[2] for x in results]
#     esnull_s = [x[3].tolist() for x in results]
#     gs = gspa_significance(es_s, esnull_s)
#     df_results = pd.DataFrame([x for x in gs],columns=['es_0','nes_0','pval_0','fdr_0'])
#     df_results['gene_set_name'] = pd.DataFrame([x[0] for x in results], index=df_results.index)
#     df_results['size0'] = pd.DataFrame([x[1] for x in results], index=df_results.index)
#     df_results['size1'] = pd.DataFrame([x[10] for x in results], index=df_results.index)

#     df_results.set_index('gene_set_name',inplace=True,drop=True)

#     es_s = [x[6] for x in results]
#     esnull_s = [x[7].tolist() for x in results]
#     gs = gspa_significance(es_s, esnull_s)
#     df_results[['es_1','nes_1','pval_1','fdr_1']] = pd.DataFrame([x for x in gs], index=df_results.index)
#     # df_results['size'] = pd.DataFrame([x[1] for x in results], index=df_results.index)

#     df_results.sort_values(by='fdr_1',ascending=True,inplace=True)

#     df_results.to_csv('results_28Jan2022/'+name+'.csv')
# # %%
