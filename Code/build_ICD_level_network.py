#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 14:51:23 2020

@author: beatriz
"""

import os
#import glob
import numpy as np
import pandas as pd
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

os.chdir("/home/beatriz/Desktop/Desktop_20200508/Network_building")

results_dir = "../ICD9_RESULTS_GREIN/"
all_files = os.listdir(results_dir)
print(len(all_files))
unwanted = {'is_finished.txt','normalized_counts','after_combat_counts','after_qrdecomp_counts','dl_general_info.txt','dl_general_info_plus_umls.txt'}
disease_dirs = set(all_files).difference(unwanted)
n_diseases = len(disease_dirs)
print(n_diseases) # 41

dm_table = pd.read_table("../../ANALYSIS/ICD9_Disease_gene_variability/all_diseases_dm_ratio.csv") # UNCOMMENT FOR DM ANALYSIS

print(disease_dirs)

# Create dictionaries with the disease information --------------------------
allgenes_dic = {}   # KEY: disease name  VALUE: list of all genes
sDEG_dic = {}       # KEY: disease name  VALUE: list of sDEGs
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
    #    dm_dic[dis] = dm_table[['genes',dis]]
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

    
# Build df will all the sDEGs / all genes with a 0 logFC - to complete expression values of a given disease
zero_sgenes_df = pd.DataFrame(dict(zip(list(sgene_union),list(np.zeros(len(sgene_union))))),index=[0])
zero_sgenes_df = zero_sgenes_df.transpose() ;  zero_sgenes_df.reset_index(inplace=True) ; 
zero_sgenes_df.rename(columns={'index':'symbol',0:'logFC'}, inplace = True)
zero_allgenes_df = pd.DataFrame(dict(zip(list(allgene_union),list(np.zeros(len(allgene_union))))),index=[0])
zero_allgenes_df = zero_allgenes_df.transpose() ;  zero_allgenes_df.reset_index(inplace=True) ; 
zero_allgenes_df.rename(columns={'index':'symbol',0:'logFC'}, inplace = True)
    
   
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
       
       df1 = dm_dic[dis1]
       expr1 = df1[df1['genes'].isin(cintersect)]
       expr1 = expr1.sort_values(by='genes', inplace=False)
       
       df2 = df_dic[dis2]
       expr2 = df2[df2['genes'].isin(cintersect)]
       expr2 = expr2.sort_values(by='genes', inplace=False)
       
       if(len(expr1['genes']) > 1):
           val,pv = spearmanr(list(expr1[dis1]), list(expr2[dis2]))
       else:
           val,pv = (np.nan,np.nan)
       
       return(val,pv)


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
            outpath = 'icd9_distances/'+ distance_func.__name__ + "_" + expr_or_var + '.csv'
            outdf.to_csv(outpath, index=False)
        return (matrix,pv_matrix,outdf)
        
        
    elif (distance_func.__name__ in valid_distance_funcs):
        print('Other function')
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
            outpath = 'icd9_distances/'+ distance_func.__name__ + "_" + expr_or_var + '.csv'
            outdf.to_csv(outpath, index=False)
        return (matrix,outdf)
    
    else:
        raise ValueError("Non valid distance_func")
        

# GENERATING THE NETWORKS WITH THEIR DISTANCES ---------------------------------------------
sDEGs_matrix, pvmatrix, outdf = build_distance_matrix(dis_list, comparable_spearman_distance)
sDEGs_matrix, pvmatrix, outdf = build_distance_matrix(dis_list, pairwise_spearman_distance)
sDEGs_matrix, pvmatrix, outdf = build_distance_matrix(dis_list, pairwise_union_spearman_distance) 


# Testing the distance functions ---------------------------------------------
dis1='371'
dis2='742'   

pairwise_spearman_distance(dis1,dis2)
pairwise_union_spearman_distance(dis1,dis2)
comparable_spearman_distance(dis1,dis2)

pairwise_spearman_distance(dis1,dis2,genes='allgenes')
pairwise_union_spearman_distance(dis1,dis2,genes='allgenes')
comparable_spearman_distance(dis1,dis2,genes='allgenes')

sDEGs_matrix, pvmatrix, outdf = build_distance_matrix(dis_list, pairwise_spearman_distance)
sDEGs_matrix, pvmatrix, outdf = build_distance_matrix(dis_list, pairwise_union_spearman_distance)
sDEGs_matrix, pvmatrix, outdf = build_distance_matrix(dis_list, comparable_spearman_distance)

  



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