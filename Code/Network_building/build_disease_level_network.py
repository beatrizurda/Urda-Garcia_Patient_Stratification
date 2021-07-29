# -*- coding: utf-8 -*-
"""
ANALYSIS OF DISEASE SIMILARITY BASED ON (s)DEGs

1- Use different metrics to obtain the pairwise distance between the diseases
   and build distance matrix
   
2- Use k-mean clutering and Silhouette k estimation to obtain disease clustering

Beatriz Urda García       
January 2020
"""
import os
#import glob
import numpy as np
import pandas as pd
import statistics
from scipy.stats import pearsonr,spearmanr
from scipy.spatial import distance
from jaccard_index.jaccard import jaccard_index
from sklearn.cluster import KMeans
from sklearn.cluster import AffinityPropagation
from sklearn import metrics
from itertools import chain
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.metrics import silhouette_samples, silhouette_score
#from sklearn_extra.cluster import KMedoids

######## CHOOSE if the analysis is done based on DEGs
analyze_by_variability = False

%matplotlib inline
from mpl_toolkits.mplot3d import Axes3D
plt.rcParams['figure.figsize'] = (16, 9)
plt.style.use('ggplot')

os.chdir("/home/beatriz/Desktop/ANALYSIS/Network_building")

results_dir = "../DL_RESULTS_GREIN/"
all_files = os.listdir(results_dir)
print(len(all_files))
unwanted = {'is_finished.txt','normalized_counts','after_combat_counts',
            'after_qrdecomp_counts','dl_general_info.txt',
            'dl_general_info_plus_umls.txt','ranked_gene_lists'}
disease_dirs = set(all_files).difference(unwanted)
n_diseases = len(disease_dirs)
print(n_diseases)
dm_table = pd.read_table("../Disease_gene_variability/all_diseases_dm_ratio.csv") # UNCOMMENT FOR DM ANALYSIS

#disease_dirs.remove('is_finished.txt','normalized_counts','after_combat_counts','after_qrdecomp_counts','dl_general_info.txt')
print(disease_dirs)

len(disease_dirs)

# Create dictionaries with the disease information --------------------------
allgenes_dic = {}   # KEY: disease name  VALUE: list of all genes
sDEG_dic = {}       # KEY: disease name  VALUE: list of sDEGs
#sDVG_dic = {}       # KEY: disease name  VALUE: list of sDVGs
df_dic = {}         # KEYS: disease      VALUE: its DE table

dis_list = list(disease_dirs) # disease list
for k in range(0,len(disease_dirs)):
    dis = dis_list[k]
    if analyze_by_variability == False:
        print("Computing distances based on DEGs")
        df = pd.read_table(results_dir+dis+'/'+dis+'_DEGs.txt') 
        filtered = df[df['adj.P.Val'] <= 0.05]
        allgenes_dic[dis] = set(df['symbol'])
        sDEG_dic[dis] = set(filtered['symbol'])
        df_dic[dis] = df
    else:
        pass
    
    
# Obtain the UNION of genes and sDEGs for ALL diseases
allgene_union = set(chain(*allgenes_dic.values())) # 20020
sgene_union = set(chain(*sDEG_dic.values())) # 17063
 

# Distance functions ---------------------------------------------------------
valid_distance_funcs = ['overlapping_distance','jaccardindex_distance',
                        'pearson_distance','spearman_distance',
                        'hamming_distance','ponderated_hamming_distance',
                        'comparable_hamming_distance',
                        'comparable_spearman_distance',
                        'pairwise_spearman_distance',
                        'pairwise_union_spearman_distance']

def get_overlap(dis_list):
    
    # Create an empty matrix
    matrix = np.zeros((n_diseases,n_diseases))
    overlapl = list()
    
    for k1 in range(0,len(dis_list)):
        for k2 in range(k1+1,len(dis_list)):
            cintersect = sDEG_dic[dis_list[k1]].intersection(sDEG_dic[dis_list[k2]])
            if(cintersect == set()):
                matrix[k1, k2] = 0
                overlapl.append(0)
            else:
                matrix[k1, k2] = len(cintersect)
                overlapl.append(len(cintersect))
    return(matrix,overlapl)
    

overlap_m, overlap_l = get_overlap(dis_list)
overlap_m.mean() 
statistics.mean(overlap_l) # DVGS: 17.351   DEGs:  345.1353
statistics.median(overlap_l) # DVGS: 2  DEGs: 0.0

if analyze_by_variability:
    outpath = '../Overlapping_number_of_sDVGs.txt'
else:
    outpath = '../Overlapping_number_of_sDEGs.txt'
overlap_df = pd.DataFrame(overlap_l, columns=['Overlaping_genes'])
overlap_df.to_csv(outpath, index=False)

def count_significant_genes(dis_list):
    
    # Create an empty matrix
    matrix = np.zeros((n_diseases,3))
    
    for k in range(0,len(dis_list)):
        cdf = df_dic[dis_list[k]]
        filtered = cdf[cdf['adj.P.Val'] <= 0.05]
        matrix[k,1] = (filtered[filtered['logFC'] >= 0]).shape[0]
        matrix[k,2] = (filtered[filtered['logFC'] < 0]).shape[0]
        matrix[k,0] = matrix[k,1] + matrix[k,2]
        
    df = pd.DataFrame(matrix, columns=['Total','Positive','Negative'])
    df.insert(0,'disease_name',dis_list)
    df = df.sort_values(by=['disease_name'])
    if analyze_by_variability:
        outpath = '../Number_over_under_sDVGs.txt'
    else:
        outpath = '../Number_over_under_sDEGs.txt'  
    df.to_csv(outpath, index=False, sep='\t')
    return(df)

table_over_under_genes = count_significant_genes(dis_list)

def jaccard_index(set1,set2):
    '''Returns the jaccard index of two sets'''
    
    intersect_len = len(set1.intersection(set2))
    union_len = len(set1.union(set2))
    if union_len == 0:
        return np.nan
    else:
        return (intersect_len/union_len)
    
def jaccardindex_distance(dis1,dis2,genes='sDEGs'):
   '''Returns the jaccard distance between two diseases based on their 
   sDEGs or allgenes
   
   jaccard_distance is 1-jaccard_index''' 
   
   if genes == 'sDEGs':
       return(1 - jaccard_index(sDEG_dic[dis1], sDEG_dic[dis2]))
   elif genes == 'allgenes':
       return(1 - jaccard_index(allgenes_dic[dis1], allgenes_dic[dis2]))
   else:
       raise ValueError("The argument genes must be 'sDEGs' or 'allgenes'")
    
def overlapping_distance(dis1,dis2,genes='sDEGs'):
   '''Number of sDEGs / all genes in common''' 
   if genes == 'sDEGs':
       return(len(sDEG_dic[dis1].intersection(sDEG_dic[dis2])))
   elif genes == 'allgenes':
       return(len(allgenes_dic[dis1].intersection(allgenes_dic[dis2])))
   else:
       raise ValueError("The argument genes must be 'sDEGs' or 'allgenes'")

def logFC_to_binary(logfc):
    '''Returns a numpy array with 3 values:
        - 1 if the value is negative
        0 if the value is 0
        1 if the value is positive
    '''
    logfc = np.where(logfc > 0, 1, logfc)
    logfc = np.where(logfc < 0, -1, logfc)
    return(logfc)
    
# Build df will all the sDEGs / all genes with a 0 logFC - to complete expression values of a given disease
zero_sgenes_df = pd.DataFrame(dict(zip(list(sgene_union),list(np.zeros(len(sgene_union))))),index=[0])
zero_sgenes_df = zero_sgenes_df.transpose() ;  zero_sgenes_df.reset_index(inplace=True) ; 
zero_sgenes_df.rename(columns={'index':'symbol',0:'logFC'}, inplace = True)
zero_allgenes_df = pd.DataFrame(dict(zip(list(allgene_union),list(np.zeros(len(allgene_union))))),index=[0])
zero_allgenes_df = zero_allgenes_df.transpose() ;  zero_allgenes_df.reset_index(inplace=True) ; 
zero_allgenes_df.rename(columns={'index':'symbol',0:'logFC'}, inplace = True)
    
def comparable_hamming_distance(dis1,dis2,genes='sDEGs'):
   '''Returns the hamming distance between 2 diseases taking into account the
   entire set of sDEGs / al genes'''
   if genes == 'sDEGs':
       df1 = df_dic[dis1]
       df1 = df1[['symbol','logFC']]
       df1 = df1[df1['symbol'].isin(sgene_union)]
       dif1 = sgene_union.difference(df1['symbol'])
       to_append = zero_sgenes_df[zero_sgenes_df['symbol'].isin(dif1)]
       df1 = df1.append(to_append, ignore_index=True) # Add the genes that are in the union of sDEGs for all diseases with a 0
       
       df2 = df_dic[dis2]
       df2 = df2[['symbol','logFC']]
       df2 = df2[df2['symbol'].isin(sgene_union)]
       dif2 = sgene_union.difference(df2['symbol'])
       to_append = zero_sgenes_df[zero_sgenes_df['symbol'].isin(dif2)]
       df2 = df2.append(to_append, ignore_index=True) # Add the genes that are in the union of sDEGs for all diseases with a 0  
       
   elif genes == 'allgenes':
       df1 = df_dic[dis1]
       df1 = df1[['symbol','logFC']]
       df1 = df1[df1['symbol'].isin(allgene_union)]
       dif1 = allgene_union.difference(df1['symbol'])
#       to_append = pd.DataFrame(dict(zip(list(dif1),list(np.zeros(len(dif1))))),index=[0])
#       to_append = to_append.transpose() ;  to_append.reset_index(inplace=True) ; 
#       to_append.rename(columns={'index':'symbol',0:'logFC'}, inplace = True)
       to_append = zero_allgenes_df[zero_allgenes_df['symbol'].isin(dif1)]
       df1 = df1.append(to_append, ignore_index=True) # Add the genes that are in the union of sDEGs for all diseases with a 0
       
       df2 = df_dic[dis2]
       df2 = df2[['symbol','logFC']]
       df2 = df2[df2['symbol'].isin(allgene_union)]
       dif2 = allgene_union.difference(df2['symbol'])
#       to_append = pd.DataFrame(dict(zip(list(dif2),list(np.zeros(len(dif2))))),index=[0])
#       to_append = to_append.transpose() ;  to_append.reset_index(inplace=True) ; 
#       to_append.rename(columns={'index':'symbol',0:'logFC'}, inplace = True)
       to_append = zero_allgenes_df[zero_allgenes_df['symbol'].isin(dif2)]
       df2 = df2.append(to_append, ignore_index=True) # Add the genes that are in the union of sDEGs for all diseases with a 0  
       
   else:
       raise ValueError("The argument genes must be 'sDEGs' or 'allgenes'")
   
   df1 = df1.sort_values(by='symbol', inplace=False)
   df2 = df2.sort_values(by='symbol', inplace=False)
   bin_vec1 = logFC_to_binary(df1['logFC'])
   bin_vec2 = logFC_to_binary(df2['logFC'])
   
   if(len(bin_vec1) > 0):
       return distance.hamming(bin_vec1, bin_vec2)
   else:
       return np.nan

def hamming_distance(dis1,dis2,genes='sDEGs'):
   ''' '''
   if genes == 'sDEGs':
       cunion = sDEG_dic[dis1].union(sDEG_dic[dis2])
   elif genes == 'allgenes':
       cunion = allgenes_dic[dis1].union(allgenes_dic[dis2])
   else:
       raise ValueError("The argument genes must be 'sDEGs' or 'allgenes'")
       
   callgenes_intersect = allgenes_dic[dis1].intersection(allgenes_dic[dis2])  
   
   df1 = df_dic[dis1]
   expr1 = df1[df1['symbol'].isin(cunion)]
   expr1 = expr1[expr1['symbol'].isin(callgenes_intersect)] # Select the union of the genes that are in both diseases
   expr1 = expr1.sort_values(by='symbol', inplace=False)
   
   df2 = df_dic[dis2]
   expr2 = df2[df2['symbol'].isin(cunion)]
   expr2 = expr2[expr2['symbol'].isin(callgenes_intersect)]
   expr2 = expr2.sort_values(by='symbol', inplace=False)
    
   bin_vec1 = logFC_to_binary(expr1['logFC'])
   bin_vec2 = logFC_to_binary(expr2['logFC'])
   
   if(len(bin_vec1) > 0):
       return distance.hamming(bin_vec1, bin_vec2)
   else:
       return np.nan
       

def ponderated_hamming_distance(dis1,dis2,genes='sDEGs'):
   ''' '''
   if genes == 'sDEGs':
       cunion = sDEG_dic[dis1].union(sDEG_dic[dis2])
   elif genes == 'allgenes':
       cunion = allgenes_dic[dis1].union(allgenes_dic[dis2])
   else:
       raise ValueError("The argument genes must be 'sDEGs' or 'allgenes'")
       
   callgenes_intersect = allgenes_dic[dis1].intersection(allgenes_dic[dis2])  
   
   df1 = df_dic[dis1]
   expr1 = df1[df1['symbol'].isin(cunion)]
   expr1 = expr1[expr1['symbol'].isin(callgenes_intersect)] # Select the union of the genes that are in both diseases
   expr1 = expr1.sort_values(by='symbol', inplace=False)
   
   df2 = df_dic[dis2]
   expr2 = df2[df2['symbol'].isin(cunion)]
   expr2 = expr2[expr2['symbol'].isin(callgenes_intersect)]
   expr2 = expr2.sort_values(by='symbol', inplace=False)
    
   bin_vec1 = logFC_to_binary(expr1['logFC'])
   bin_vec2 = logFC_to_binary(expr2['logFC'])
   
   if(len(bin_vec1) > 0):
       return (distance.hamming(bin_vec1, bin_vec2) / len(bin_vec1))
   else:
       return np.nan

   
def pearson_distance(dis1,dis2,genes='sDEGs'):
   '''
   For each disease pair, I want to keep the genes that are un the union of the 
   sDEGs but also in the intersection of the union of all genes.
   ''' 
   if genes == 'sDEGs':
       cunion = sDEG_dic[dis1].union(sDEG_dic[dis2])  
   
   elif genes == 'allgenes':
       cunion = allgenes_dic[dis1].union(allgenes_dic[dis2])
   
   else:
       raise ValueError("The argument genes must be 'sDEGs' or 'allgenes'")
   
   callgenes_intersect = allgenes_dic[dis1].intersection(allgenes_dic[dis2])    
   
   df1 = df_dic[dis1]
   expr1 = df1[df1['symbol'].isin(cunion)]
   expr1 = expr1[expr1['symbol'].isin(callgenes_intersect)] # Select the union of the genes that are in both diseases
   expr1 = expr1.sort_values(by='symbol', inplace=False)
   
   df2 = df_dic[dis2]
   expr2 = df2[df2['symbol'].isin(cunion)]
   expr2 = expr2[expr2['symbol'].isin(callgenes_intersect)]
   expr2 = expr2.sort_values(by='symbol', inplace=False)
   
   if(len(expr1['symbol']) > 1):
       val,pv = pearsonr(list(expr1['logFC']), list(expr2['logFC']))
   else:
       val,pv = (np.nan,np.nan)
   
   return(val,pv)
   
def spearman_distance(dis1,dis2,genes='sDEGs'):
   ''' ''' 
   if genes == 'sDEGs':
       cunion = sDEG_dic[dis1].union(sDEG_dic[dis2])  
   
   elif genes == 'allgenes':
       cunion = allgenes_dic[dis1].union(allgenes_dic[dis2])
   
   else:
       raise ValueError("The argument genes must be 'sDEGs' or 'allgenes'")
     
        
   callgenes_intersect = allgenes_dic[dis1].intersection(allgenes_dic[dis2])
       
   df1 = df_dic[dis1]
   expr1 = df1[df1['symbol'].isin(cunion)]
   expr1 = expr1[expr1['symbol'].isin(callgenes_intersect)]
   expr1 = expr1.sort_values(by='symbol', inplace=False)
   
   df2 = df_dic[dis2]
   expr2 = df2[df2['symbol'].isin(cunion)]
   expr2 = expr2[expr2['symbol'].isin(callgenes_intersect)] # Select the union of the genes that are in both diseases
   expr2 = expr2.sort_values(by='symbol', inplace=False)
   
   if(len(expr1['symbol']) > 1):
       val,pv = spearmanr(list(expr1['logFC']), list(expr2['logFC']))
   else:
       val,pv = (np.nan,np.nan)
   
   return(val,pv)
   
   
def comparable_spearman_distance(dis1,dis2,genes='sDEGs'):
   '''Returns the spearman distance between 2 diseases taking into account the
   entire set of sDEGs / al genes; that is the union of the sDEGs and all genes 
   for all diseases, putting zeros if needed'''
   if genes == 'sDEGs':
       df1 = df_dic[dis1]
       df1 = df1[['symbol','logFC']]
       df1 = df1[df1['symbol'].isin(sgene_union)]
       dif1 = sgene_union.difference(df1['symbol'])
       to_append = zero_sgenes_df[zero_sgenes_df['symbol'].isin(dif1)]
       df1 = df1.append(to_append, ignore_index=True) # Add the genes that are in the union of sDEGs for all diseases with a 0
       
       df2 = df_dic[dis2]
       df2 = df2[['symbol','logFC']]
       df2 = df2[df2['symbol'].isin(sgene_union)]
       dif2 = sgene_union.difference(df2['symbol'])
       to_append = zero_sgenes_df[zero_sgenes_df['symbol'].isin(dif2)]
       df2 = df2.append(to_append, ignore_index=True) # Add the genes that are in the union of sDEGs for all diseases with a 0  
       
   elif genes == 'allgenes':
       df1 = df_dic[dis1]
       df1 = df1[['symbol','logFC']]
       df1 = df1[df1['symbol'].isin(allgene_union)]
       dif1 = allgene_union.difference(df1['symbol'])
       to_append = zero_allgenes_df[zero_allgenes_df['symbol'].isin(dif1)]
       df1 = df1.append(to_append, ignore_index=True) # Add the genes that are in the union of sDEGs for all diseases with a 0
       
       df2 = df_dic[dis2]
       df2 = df2[['symbol','logFC']]
       df2 = df2[df2['symbol'].isin(allgene_union)]
       dif2 = allgene_union.difference(df2['symbol'])
       to_append = zero_allgenes_df[zero_allgenes_df['symbol'].isin(dif2)]
       df2 = df2.append(to_append, ignore_index=True) # Add the genes that are in the union of sDEGs for all diseases with a 0  
       
   else:
       raise ValueError("The argument genes must be 'sDEGs' or 'allgenes'")
   
   df1 = df1.sort_values(by='symbol', inplace=False)
   df2 = df2.sort_values(by='symbol', inplace=False)
   
   if(len(df1['symbol']) > 1):
       val,pv = spearmanr(list(df1['logFC']), list(df2['logFC']))
   else:
       val,pv = (np.nan,np.nan)
   
   return(val,pv)


def pairwise_spearman_distance(dis1,dis2,genes='sDEGs',comparison='intersection'):
   ''' ''' 
#   dis1 = 'Asthma'
#   dis2 = 'Asthma'
#   dis2 = 'Schizophrenia'
#   genes='sDEGs'
#   genes='DM'
   
   if genes == 'sDEGs':
       cintersect = sDEG_dic[dis1].intersection(sDEG_dic[dis2]) 
   
   elif ((genes == 'allgenes') or (genes == 'DM')):
       cintersect = allgenes_dic[dis1].intersection(allgenes_dic[dis2])
   
   else:
       raise ValueError("The argument genes must be 'sDEGs', 'allgenes' or 'DM'")
     
       
   if ((genes == 'sDEGs') or (genes == 'allgenes')):
       
       df1 = df_dic[dis1]
       expr1 = df1[df1['symbol'].isin(cintersect)]
       expr1 = expr1.sort_values(by='symbol', inplace=False)
       
       df2 = df_dic[dis2]
       expr2 = df2[df2['symbol'].isin(cintersect)]
       expr2 = expr2.sort_values(by='symbol', inplace=False)
       
       if(len(expr1['symbol']) > 1):
           val,pv = spearmanr(list(expr1['logFC']), list(expr2['logFC']))
       else:
           val,pv = (np.nan,np.nan)
       
       return(val,pv)
       
   elif genes == 'DM':
       
       pass


def pairwise_union_spearman_distance(dis1,dis2,genes='sDEGs'):
   ''' ''' 
#   dis1 = 'Asthma'
#   dis2 = 'Asthma'
#   dis2 = 'Schizophrenia'
#   genes='sDEGs'
#   genes='DM'
   
   if genes == 'sDEGs':
       cunion = sDEG_dic[dis1].union(sDEG_dic[dis2]) 
   
   elif ((genes == 'allgenes') or (genes == 'DM')):
       cunion = allgenes_dic[dis1].union(allgenes_dic[dis2])
   
   else:
       raise ValueError("The argument genes must be 'sDEGs', 'allgenes' or 'DM'")
     
       
   if ((genes == 'sDEGs') or (genes == 'allgenes')):
       
       df1 = df_dic[dis1] ; df1 = df1[['symbol','logFC']]
       expr1 = df1[df1['symbol'].isin(cunion)]
       # Add the genes that are in the union, with their logFC value or 0 if they are lowly expressed
       dif1 = cunion.difference(expr1['symbol'])
       to_append = zero_allgenes_df[zero_allgenes_df['symbol'].isin(dif1)]
       expr1 = expr1.append(to_append, ignore_index=True)
       expr1 = expr1.sort_values(by='symbol', inplace=False)
       
       df2 = df_dic[dis2] ; df2 = df2[['symbol','logFC']]
       expr2 = df2[df2['symbol'].isin(cunion)]
       # Add the genes that are in the union, with their logFC value or 0 if they are lowly expressed
       dif2 = cunion.difference(expr2['symbol'])
       to_append = zero_allgenes_df[zero_allgenes_df['symbol'].isin(dif2)]
       expr2 = expr2.append(to_append, ignore_index=True)
       expr2 = expr2.sort_values(by='symbol', inplace=False)
       
       if(len(expr1['symbol']) > 1):
           val,pv = spearmanr(list(expr1['logFC']), list(expr2['logFC']))
       else:
           val,pv = (np.nan,np.nan)
       
       return(val,pv)
       
   elif genes == 'DM':
       pass
            
   
# Create the distance matrix ------------------------------------------------- 
def build_distance_matrix(dis_list,distance_func, genes='sDEGs', sort=True, save_results=True):
    ''''''
    # Create an empty matrix
    matrix = np.zeros((n_diseases,n_diseases))
    
    if ((distance_func.__name__ == 'pearson_distance') or (distance_func.__name__ == 'spearman_distance') or (distance_func.__name__ == 'comparable_spearman_distance') or (distance_func.__name__ == 'pairwise_union_spearman_distance') or (distance_func.__name__ == 'pairwise_spearman_distance')):
        outdf = pd.DataFrame(columns=['Dis1','Dis2','Distance','pvalue']) #output df
        pv_matrix = np.zeros((n_diseases,n_diseases))
       
        for k1 in range(0,len(dis_list)):
            for k2 in range(k1+1,len(dis_list)):
                try:
#                    print(dis_list[k1] + '  ' + dis_list[k2])
                    cdis,cpv = distance_func(dis_list[k1],dis_list[k2],genes=genes)
    #                print(dis_list[k1]+" "+dis_list[k2]+" --> "+str('%.4f' %cdis)+"  pv= "+('%.4f' %cpv))
                    matrix[k1, k2] = cdis
                    pv_matrix[k1, k2] = cpv
                    outdf = outdf.append({'Dis1':dis_list[k1],'Dis2':dis_list[k2],'Distance':cdis,'pvalue':cpv},ignore_index=True)
                except:
                    print(dis_list[k1] + ' ' + dis_list[k2] + ' error')
                
        if sort:
            outdf = outdf.sort_values(by=['Distance','pvalue'], ascending=False, inplace=False)
        if save_results:
            if analyze_by_variability == True:
                expr_or_var = 'sDVGs'
            else:
                expr_or_var = genes
            outpath = 'distances/'+ distance_func.__name__ + "_" + expr_or_var + '.csv'
            outdf.to_csv(outpath, index=False)
        return (matrix,pv_matrix,outdf)
        
        
    elif (distance_func.__name__ in valid_distance_funcs):
        outdf = pd.DataFrame(columns=['Dis1','Dis2','Distance']) #output df
        for k1 in range(0,len(dis_list)):
            for k2 in range(k1+1,len(dis_list)):
                cdis = distance_func(dis_list[k1],dis_list[k2],genes=genes)
#                print(dis_list[k1]+" "+dis_list[k2]+" --> "+str(cdis))
                matrix[k1, k2] = cdis
                outdf = outdf.append({'Dis1':dis_list[k1],'Dis2':dis_list[k2],'Distance':cdis},ignore_index=True)
                
        if sort:
            outdf = outdf.sort_values(by='Distance', ascending=False, inplace=False, na_position='last')
        if save_results:
            if analyze_by_variability == True:
                expr_or_var = 'sDVGs'
            else:
                expr_or_var = genes
            outpath = 'distances/'+ distance_func.__name__ + "_" + expr_or_var + '.csv'
            outdf.to_csv(outpath, index=False)
        return (matrix,outdf)
    
    else:
        raise ValueError("Non valid distance_func")
 

# GENERATING THE NETWORKS WITH THEIR DISTANCES ---------------------------------------------
sDEGs_matrix, pvmatrix, outdf = build_distance_matrix(dis_list, comparable_spearman_distance)
sDEGs_matrix, pvmatrix, outdf = build_distance_matrix(dis_list, pairwise_spearman_distance)
sDEGs_matrix, pvmatrix, outdf = build_distance_matrix(dis_list, pairwise_union_spearman_distance)      
        
# Testing the distance functions ---------------------------------------------
dis1='LungCancer'
dis2='Asthma'      
print(overlapping_distance(dis1,dis1))
print(pearson_distance(dis1,dis1))
print(spearman_distance(dis1,dis1))
print(spearman_distance(dis1,dis2))
print(spearman_distance(dis1,'Asthma', genes='sDEGs'))
print(spearman_distance(dis1,'Asthma',genes='allgenes'))
print(pearson_distance(dis1,'Asthma', genes='sDEGs'))   ##### I DON't get right the all genes.
print(pearson_distance(dis1,'Asthma',genes='allgenes'))
print(jaccardindex_distance(dis1,dis1))
print(jaccardindex_distance(dis1,'Asthma'))
print(pairwise_spearman_distance(dis1,'Asthma'))
hamming_distance(dis1,dis2,genes='allgenes')
ponderated_hamming_distance(dis1,dis2,genes='allgenes')
comparable_hamming_distance(dis1,dis2,genes='allgenes')

comparable_spearman_distance(dis1,dis2,genes='allgenes')
comparable_spearman_distance(dis1,dis2)
pairwise_union_spearman_distance(dis1,dis2)
pairwise_union_spearman_distance(dis1,dis2,genes='allgenes')

sDEGs_matrix, pvmatrix, outdf = build_distance_matrix([dis1,'Asthma'], pearson_distance)
sDEGs_matrix, pvmatrix, outdf = build_distance_matrix([dis1,'Asthma'], pearson_distance, genes='allgenes')

# Creating matrices with sDEGs -----------------------------------------------
sDEGs_matrix, outdf = build_distance_matrix(dis_list, overlapping_distance)
jacc_sDEGs_matrix, outdf = build_distance_matrix(dis_list, jaccardindex_distance)
sDEGs_matrix, outdf = build_distance_matrix(dis_list, hamming_distance)
sDEGs_matrix, outdf = build_distance_matrix(dis_list, ponderated_hamming_distance)
sDEGs_matrix, pvmatrix, outdf = build_distance_matrix(dis_list, pearson_distance)
sDEGs_matrix, pvmatrix, outdf = build_distance_matrix(dis_list, spearman_distance)
sDEGs_matrix, outdf = build_distance_matrix(dis_list, comparable_hamming_distance)

sDEGs_matrix, pvmatrix, outdf = build_distance_matrix(dis_list, comparable_spearman_distance)
sDEGs_matrix, pvmatrix, outdf = build_distance_matrix(dis_list, pairwise_spearman_distance)
sDEGs_matrix, pvmatrix, outdf = build_distance_matrix(dis_list, pairwise_union_spearman_distance)
sDEGs_matrix, pvmatrix, outdf = build_distance_matrix(['Ulcer','Ischemia'], pairwise_union_spearman_distance)
pairwise_union_spearman_distance('Ulcer','Ischemia',genes='allgenes')
sDEGs_matrix, pvmatrix, outdf = build_distance_matrix(dis_list, pairwise_spearman)

# Creating matrices with all genes -------------------------------------------
sDEGs_matrix, outdf = build_distance_matrix(dis_list, overlapping_distance, genes='allgenes')
sDEGs_matrix, outdf = build_distance_matrix(dis_list, hamming_distance, genes='allgenes')
sDEGs_matrix, outdf = build_distance_matrix(dis_list,  ponderated_hamming_distance, genes='allgenes')
sDEGs_matrix, pvmatrix, outdf = build_distance_matrix(dis_list, pearson_distance, genes='allgenes')
sDEGs_matrix, pvmatrix, outdf = build_distance_matrix(dis_list, spearman_distance, genes='allgenes')
all_matrix, outdf = build_distance_matrix(dis_list, comparable_hamming_distance, genes='allgenes')

sDEGs_matrix, pvmatrix, outdf = build_distance_matrix(dis_list, comparable_spearman_distance, genes='allgenes')

### Exporting matrices to plot them in R
jacdf = pd.DataFrame(data = jacc_sDEGs_matrix, index=dis_list, columns=dis_list)
ham_sde = pd.DataFrame(data = sDEGs_matrix, index=dis_list, columns=dis_list)
ham_allg = pd.DataFrame(data = all_matrix, index=dis_list, columns=dis_list)

jacdf = jacdf.fillna(0)
jacdf.to_csv('diseases_jaccard_distance_sDEGs.csv', index=True)
ham_sde.to_csv('diseases_hamming_distance_sDEGs.csv', index=True)
ham_allg.to_csv('diseases_hamming_distance_all_genes.csv', index=True)


# Creating one df and file with everything - all metrics for each disease pair -----------------
distance_funcs = [pearson_distance, spearman_distance, 
                  hamming_distance, ponderated_hamming_distance,
                  overlapping_distance, jaccardindex_distance,
                  comparable_hamming_distance,
                  comparable_spearman_distance]
distance_funcs = [spearman_distance, comparable_hamming_distance,comparable_spearman_distance]
genes_opt = ['sDEGs','allgenes']
entered = False
for k in distance_funcs:
    cname = k.__name__.replace("_distance","")
    for j in genes_opt:
        if ((k.__name__ == 'pearson_distance') or (k.__name__ == 'spearman_distance')or (k.__name__ == 'comparable_spearman_distance')):
            first, second, coutdf = build_distance_matrix(dis_list, k, j,sort=False, save_results=False)
            if entered == False:
                coutdf.columns.values[2] = cname+"_"+j
                finaldf = coutdf
                entered = True
                print(coutdf.shape)
            else:
                print(coutdf.shape)
                finaldf[str(cname+"_"+j)] = coutdf['Distance'].values
                finaldf[str(cname+'pvalue')] = coutdf['pvalue'].values
#                coutdf = coutdf.drop(labels=['Dis1','Dis2'])
#                finaldf = pd.concat([finaldf,coutdf], axis=1)
        else:
            fist, coutdf = build_distance_matrix(dis_list, k, j,sort=False, save_results=False)
            finaldf[str(cname+"_"+j)] = coutdf['Distance'].values
        
finaldf = finaldf.sort_values(by='spearman_sDEGs', ascending=False, inplace=False, na_position='last')
finaldf.to_csv('disease_distances_all_metrics.csv', index=False)  

finaldf = finaldf.sort_values(by='spearman_sDEGs', ascending=False, inplace=False, na_position='last')
finaldf.to_csv('disease_distances_imp_metrics.csv', index=False)

# Computing the correlations between the different metrics -------------------
spearmanr(list(finaldf['pearson_sDEGs']), list(finaldf['spearman_sDEGs'])) # 0.9630737525486307, pvalue=0.0
spearmanr(list(finaldf['pearson_sDEGs']), list(finaldf['overlapping_sDEGs'])) # 0.15983541818316466, pvalue=4.84335559768268e-09
spearmanr(list(finaldf['pearson_sDEGs']), list(finaldf['jaccardindex_sDEGs'])) # correlation=nan, pvalue=nan

spearmanr(list(finaldf['spearman_sDEGs']), list(finaldf['overlapping_sDEGs'])) # 0.14201802596035865, pvalue=2.0698883319417267e-07
spearmanr(list(finaldf['spearman_sDEGs']), list(finaldf['jaccardindex_sDEGs'])) 

spearmanr(list(finaldf['overlapping_sDEGs']), list(finaldf['jaccardindex_sDEGs']))

##### K-MEANS CLUSTERING OF THE DISEASES --------------------------------------------------

###### GENERATE TABLES FOR DISEASE CLUSTERING AND DISEASE SIMILARITY HEATMAP #---START---------##### 
# Obtain the UNION of genes and sDEGs for ALL diseases
allgene_union = set(chain(*allgenes_dic.values())) # 20020
sgene_union = set(chain(*sDEG_dic.values())) # 17063

# Build the dataframe to use
#alldf = pd.DataFrame(index=list(allgene_union), columns=dis_list) # all genes df -- 20020
#sdf = pd.DataFrame(index=list(sgene_union), columns=dis_list)   # sDEGs df -- 17063

alldf = pd.DataFrame(index=list(allgene_union)) # all genes df -- 20020
sdf = pd.DataFrame(index=list(sgene_union))   # sDEGs df -- 17063

alldf.shape # 19992
sdf.shape # 16865

alldf['symbol'] = list(allgene_union)
sdf['symbol'] = list(sgene_union)


# alldf.info()
#sdf.info()

# Fill the dataframes -- BEFORE: in the sdf only putting sDEGs
#for dis in dis_list:
#    cdf = df_dic[dis].copy()
#    cdf = cdf[['symbol','logFC']]
#    cdf.columns.values[1] = dis
#    csdf = cdf[cdf['symbol'].isin(sDEG_dic[dis])]
#
#    alldf = alldf.merge(cdf,how='left',on='symbol')
#    sdf = sdf.merge(csdf,how='left',on='symbol')
    
# Fill the dataframes -- NOW: putting the sDEGs that are in the union of sDEGs
for dis in dis_list:
    cdf = df_dic[dis].copy()
    cdf = cdf[['symbol','logFC']]
    cdf.columns.values[1] = dis
    csdf = cdf[cdf['symbol'].isin(sgene_union)]

    alldf = alldf.merge(cdf,how='left',on='symbol')
    sdf = sdf.merge(csdf,how='left',on='symbol')
    

#### <TO SAVE the newly generated tables> ###########################
alldf = alldf.fillna(0)
sdf = sdf.fillna(0)

alldf = alldf.set_index('symbol')
sdf = sdf.set_index('symbol')

alldf = alldf.transpose()
sdf = sdf.transpose()
alldf.to_csv('all_genes_all_diseases_logFC_table.csv', index=True)
sdf.to_csv('union_sDEGs_all_diseases_logFC_table.csv', index=True)

#alldf.rename(columns=alldf.iloc[0,], inplace=True)
#sdf.rename(columns=sdf.iloc[0,], inplace=True)

#### </> #############################################################
###### GENERATE TABLES FOR DISEASE CLUSTERING AND DISEASE SIMILARITY HEATMAP #---END---------##### 


alldf.to_csv('all_genes_all_diseases_logFC_table.csv', index=False)
sdf.to_csv('union_sDEGs_all_diseases_logFC_table.csv', index=False)
    
# <Checking that the tables have been generated correctly>   
hey = sdf[np.isnan(sdf['Autism']) == False]
hey[['symbol','Autism']]

sDEG_dic['Autism']
# </> All good :)

# Changing the symbol position from column to index
alldf = alldf.set_index('symbol')
sdf = sdf.set_index('symbol')

#disease_info = pd.read_csv("../Patients_information_table.txt", sep="\t")
disease_info = pd.read_csv("../DL_RESULTS_GREIN/dl_general_info_plus_umls.txt", sep="\t")


#umls_info <- read.csv2("Patients_information_table.txt",header=T, sep="\t",stringsAsFactors = F, row.names = NULL)

#### Perform K-MEANS algorithm

# Transpose matrix
X = alldf.transpose()
#X = sdf.transpose()
X = X.fillna(0)

# Normalize it
xx = X.values
min_max_scaler = preprocessing.MinMaxScaler()
x_scaled = min_max_scaler.fit_transform(xx)
Xn = pd.DataFrame(x_scaled,index=X.index , columns=X.columns)

# Set seed
np.random.seed(1)

# Sort the metadata information in the same order
dismeta = pd.read_csv("../DL_RESULTS_GREIN/dl_general_info_plus_umls.txt", sep="\t")
dismeta = dismeta.set_index('disease_name')
dismeta = dismeta.reindex(index=X.index)
dismeta = dismeta.reset_index()
dismeta.columns.values[0] = 'disease_name' 

# DEFINE COLOURS
# Define colors for the disease category
dis_cats = np.unique(dismeta['disease_cat'])
n_dis_cats = len(dis_cats)
discat_cols = ['blueviolet','red','springgreen1','green4','magenta','blue','brown','orange','pink2','grey']
#discat_cols = ['blue','green','red','cyan','magenta','yellow','black','orange','pink','grey']
discat_colors_dic = dict(zip(dis_cats,discat_cols))
discat_cols = list()
for k in range(0,n_diseases):
    cdiscat = dismeta['disease_cat'].values[k]
    discat_cols.append(discat_colors_dic[cdiscat])

dismeta['dis_cat_colors'] = discat_cols

# Define colors for the tissue category
tis_cats = np.unique(dismeta['tissue_cat'])
n_tis_cats = len(tis_cats)
tiscat_cols = ['black','grey','fuchsia','lightcoral','brown','red','coral','peru',
               'royalblue','gold','darkkhaki','olivedrab','greenyellow','springgreen',
               'aquamarine','mediumturquoise','darkcyan','deepskyblue','blue','blueviolet','hotpink']

tiscat_colors_dic = dict(zip(tis_cats,tiscat_cols))
tiscat_cols = list()
for k in range(0,n_diseases):
    ctiscat = dismeta['tissue_cat'].values[k]
    tiscat_cols.append(tiscat_colors_dic[ctiscat])

dismeta['tis_cat_colors'] = tiscat_cols   

# Use Silhouette to obtain the best k
silhouette = pd.DataFrame(columns=['N_clusters','Silhouette'])
for n_clusters in range(2,52):
   clusterer = KMeans(n_clusters=n_clusters, random_state=10)
   cluster_labels = clusterer.fit_predict(Xn)
   silhouette_avg = silhouette_score(Xn, cluster_labels)
   silhouette = silhouette.append({'N_clusters':n_clusters,'Silhouette':silhouette_avg},ignore_index=True)
#   silhouette.append([n_clusters,silhouette_avg])
   print("For n_clusters =", n_clusters,
          "The average silhouette_score is :", silhouette_avg)

plt.plot(silhouette['N_clusters'],silhouette['Silhouette'])
plt.show()

allbestk_index = silhouette['Silhouette'].idxmax()
allbestk_info = silhouette.iloc[allbestk_index] # 3 clusters with a Silhouette = 0.203748
allbestk = int(allbestk_info['N_clusters'])

# Compute k-means
Xn = np.array(Xn)
Xn.shape # (52, 17063)
kmeans = KMeans(n_clusters=allbestk).fit(Xn)
centroids = kmeans.cluster_centers_
print(centroids)

labels = kmeans.predict(Xn)
C = kmeans.cluster_centers_
colors = dismeta['dis_cat_colors'].values    # To colour by disease category
#colors = dismeta['tis_cat_colors'].values    # To color by tissue category
#assign = ['blue','green','red','black','magenta','orange']
assign = ['blue','green','red','cyan','magenta','yellow','black','orange','pink','grey'] 
assign = assign[0:allbestk]  
fig = plt.figure()
ax = Axes3D(fig)

# In 3D
ax.scatter(Xn[:, 0], Xn[:, 1], Xn[:, 2], c=colors,s=60)
ax.scatter(C[:, 0], C[:, 1], C[:, 2], marker='*', c=assign, s=1000)

# In 2D
plt.scatter(Xn[:, 0], Xn[:, 1], c=colors)
plt.show()

##### K-MEDOIDS CLUSTERING (PAM) --------------------------------------------
from sklearn_extra.cluster import KMedoids
kmedoids = KMedoids(n_clusters=2, random_state=0).fit(Xn)
kmedoids.labels_

from pyclustering.cluster.kmedoids import kmedoids
from pyclustering.cluster import cluster_visualizer

# Set random initial medoids.
initial_medoids = [1, 500]
# Create instance of K-Medoids algorithm.
kmedoids_instance = kmedoids(sample, initial_medoids)
# Run cluster analysis and obtain results.
kmedoids_instance.process()
clusters = kmedoids_instance.get_clusters()
# Show allocated clusters.
print(clusters)
# Display clusters.
visualizer = cluster_visualizer()
visualizer.append_clusters(clusters, sample)
visualizer.show()



##### Affinity propagation clustering ----------------------------------------
sDEGs_matrix, pvmatrix, outdf = build_distance_matrix(dis_list, pearson_distance, genes='allgenes')
X = sDEGs_matrix
#labels_true = 
af = AffinityPropagation(preference=-50).fit(X)
cluster_centers_indices = af.cluster_centers_indices_
labels = af.labels_
n_clusters_ = len(cluster_centers_indices)

print('Estimated number of clusters: %d' % n_clusters_)
print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels))
print("Completeness: %0.3f" % metrics.completeness_score(labels_true, labels))
print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels))
print("Adjusted Rand Index: %0.3f"
      % metrics.adjusted_rand_score(labels_true, labels))
print("Adjusted Mutual Information: %0.3f"
      % metrics.adjusted_mutual_info_score(labels_true, labels))

plt.close('all')
plt.figure(1)
plt.clf()
from itertools import cycle

colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
for k, col in zip(range(n_clusters_), colors):
    class_members = labels == k
    cluster_center = X[cluster_centers_indices[k]]
    plt.plot(X[class_members, 0], X[class_members, 1], col + '.')
    plt.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=14)
    for x in X[class_members]:
        plt.plot([cluster_center[0], x[0]], [cluster_center[1], x[1]], col)

plt.title('Estimated number of clusters: %d' % n_clusters_)
plt.show()

from sklearn.datasets import make_blobs
centers = [[1, 1], [-1, -1], [1, -1]]
X, labels_true = make_blobs(n_samples=300, centers=centers, cluster_std=0.5,
                            random_state=0)

#print("Silhouette Coefficient: %0.3f"
#      % metrics.silhouette_score(X, labels, metric='sqeuclidean'))




