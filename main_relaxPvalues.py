import numpy as np
from scipy.stats import norm
from scipy import integrate
from padjust import padjust
from computeK2 import computeK2
from scipy.stats import norm


def relaxPvalues(p_values,SubregionsIndices,Alpha,globalAlpha,method,q):
    '''
    This code implements the 2-step strategy from Meskaldji et al. 2015, 
    'Improved Statistical Evaluation of Group Differences in Connectomes by Screening-filtering Strategy with Application to Study Maturation of Brain 
    Connections between Childhood and Adolescence', NeuroImage, vol. 108 pp. 251-264 

    Translated to Python from Serafeim Loukas, serafeim.loukas@epfl.ch, Jul 2019.

    INPUT arguments:
        - p_values: orignial p_values (correspondingedges/nodes or any single unit) , a list of len(p_values) == len(SubregionsIndices)
        - SubregionsIndices: vector with p_value assignment. Each p_value belongs to a subset, a list of len(p_values) == len(SubregionsIndices)
        - Alpha: subset screening threshold (default 0.05). The screening threshold for subset p_values is Alpha in the Soft Thresholding Screening and Filtering (STSF) case.  The screening threshold in the Hard Thresholding SF (HTSF) is determined by the multiplicity correction used in the screening step. 
        - globalAlpha: the global type I error rate at the filtering level (default 0.05)
        - method: multiplicity correction used in the screening/filtering step ('fdr' or 'BH' for Bonferroni). Other procedures could be applied on the relaxed p_values. 
        - q: display report or not (Default true)

    RETURNS:
        - subsetPavlues: p values of subnets in numpy.ndarray format.
        - relaxedpvalues: 1x2 cell with relaxed p_values after applying the HTSF and STSF strategies in in numpy.ndarray format (result uses padjust(relaxedpvalues,method)<globalAlpha).
        - subsetScores: summary measure of subnets (uses mean, change if needed) in numpy.ndarray format.


    Example:

    from padjust import padjust
    from computeK2 import computeK2
    from main_relaxPvalues import relaxPvalues
    
    SubregionsIndices = [1,1,1,2,2,2]
    p_values = [0.001, 0.0005, 0.0001, 0.2, 0.4, 0.2]
    Alpha = 0.05
    globalAlpha = 0.05
    method = 'fdr'
    q = 1

    subsetPvalues,relaxedpvalues,subsetScores = relaxPvalues(p_values,SubregionsIndices,Alpha,globalAlpha,method,q)

    '''

    M = len(p_values) # number of p_values
    sm = len(np.unique(SubregionsIndices)) # number of subsets
    M = float(M)
    si = M / sm
    
    if not SubregionsIndices:
        subsets = np.unique(SubregionsIndices)
        if len(subsets) < 2:
            raise Exception('less than 2 subsets found')
    
    if not SubregionsIndices and (len(SubregionsIndices) != M):
        raise Exception('Number of p_values not equal to SubregionsIndices length. check data')
    
    # Subset scores and subset p_values
    scores = norm.ppf(1. - np.array(p_values))
    subsetScores = np.zeros((sm,1)) 
    subsetSizes = np.zeros((sm,1))
    for i in range(1,sm+1):
        subsetSizes[i-1,:] = np.sum(np.array(SubregionsIndices) == i)
        subsetScores[i-1,:] = (subsetSizes[i-1]**0.5) * np.mean(scores[np.array(SubregionsIndices) == i])
    subsetPvalues = 1. - norm.cdf(subsetScores)
    
    # Screening step
    R = np.sum(padjust(subsetPvalues,method) <= Alpha)
    
    if q:
        print('{} positive subsets ( {} correction)\n'.format(R,method))
        #if R > 0: print(np.sum(padjust(subsetPavlues,method) <= Alpha))
    
    ui = list()
    ui.append(norm.ppf(1.- Alpha/(sm))) # Bonferroni
    if 'fdr' in method:
        if R < 1: R = 1 # will be at least 1
        ui[0]= norm.ppf(1. - R * Alpha / (sm)) # max(1,R) will be at least 1
    ui.append(norm.ppf(1. - Alpha)) # uncorrected 
        
    
    # compute relaxation coefficient and filtering step 
    Correct = ['HTSF','STSF']
    relaxedpvalues = []
    for i in range(len(Correct)):
        relax = computeK2(globalAlpha, sm, si,ui[i])
        VectorOfSignificantSubsets= (subsetPvalues <= (1. - norm.cdf(ui[i])))
        Filter = np.zeros((int(M),1))
        for j in range(1,sm+1):
            Filter[np.array(SubregionsIndices)==j,0] = VectorOfSignificantSubsets[j-1]
        # all with nonsig subsets = 1, the others are reduced
        relaxedpvalues.append((np.array(p_values)/relax).reshape(-1,1) * Filter + (1. - Filter))
        
        temp = padjust(relaxedpvalues[i],method) <= globalAlpha
        R2 = np.sum(temp)
        if q: print('{} significant subsets and {} significant atoms (edges/nodes or any single unit) with {} ({} correction)\n'.format(np.sum(VectorOfSignificantSubsets),R2,Correct[i],method))
    
    return subsetPvalues,relaxedpvalues,subsetScores
    