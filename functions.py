
#%%
import os

import numpy as np
import pandas as pd

import pickle
import random as random

from scipy.spatial import distance

import matplotlib.pyplot as plt
import sys, logging

from math import ceil
import argparse



#%%


gene_list_input = None
gene_list_innetwork = None
d_full = None
gene_list_innetwork_indices = None
d_big = None
# 
gene_list_outnetwork_indices = None


def gmt2dict(gmt):
    """
    Load gmt file as dictionary
    """
    with open(gmt) as genesets:
        genesets_dict = { line.strip().split("\t")[0]: line.strip().split("\t")[2:]
                                    for line in genesets.readlines()}
    return genesets_dict

# 

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

# 

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
    # 
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

    # 
    Nhint = tag_indicator.sum(axis=axis, keepdims=True)
    sum_correl_tag = np.sum(correl_vector*tag_indicator, axis=axis, keepdims=True)
    # compute ES score.
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



    # 
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


#