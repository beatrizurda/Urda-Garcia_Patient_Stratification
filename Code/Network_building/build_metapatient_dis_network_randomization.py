#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 10:17:26 2021

@author: beatriz
"""

import os
#import glob
import sys
import numpy as np
import pandas as pd
from scipy.stats import pearsonr,spearmanr
from itertools import chain
from statsmodels.stats.multitest import fdrcorrection



os.chdir("Randomization/Network_building")
crandomization = sys.argv[1] # Elemento 0: nombre del script, Elemento 1: número de la randomización actual
print(crandomization)

# OPEN DISEASE FILES
results_dir = "../../DL_RESULTS_GREIN/"
all_files = os.listdir(results_dir)
print(len(all_files))
unwanted = {'is_finished.txt','normalized_counts','after_combat_counts',
            'after_qrdecomp_counts','dl_general_info.txt',
            'dl_general_info_plus_umls.txt','ranked_gene_lists'}
disease_dirs = set(all_files).difference(unwanted)
n_diseases = len(disease_dirs)
print(n_diseases)
#dm_table = pd.read_table("../Disease_gene_variability/all_diseases_dm_ratio.csv") # UNCOMMENT FOR DM ANALYSIS
print(disease_dirs)

# APPEND METAPATIENTS FILES
#m_results_dir = "../../Metapatients/with_entire_count_matrix/Randomization/DEA_results/"   # en mi ordenador
m_results_dir = "../DEA_results/"   # en el mare
m_all_files = os.listdir(m_results_dir)
print(len(m_all_files))
m_unwanted = {'is_finished.txt','normalized_counts','after_combat_counts',
            'after_qrdecomp_counts','dl_general_info.txt',
            'dl_general_info_plus_umls.txt','ranked_gene_lists'}
m_disease_dirs = set(m_all_files).difference(m_unwanted)
m_n_diseases = len(m_disease_dirs)
print(m_n_diseases)
print(m_disease_dirs)


# Create dictionaries with the disease information --------------------------
allgenes_dic = {}   # KEY: disease name  VALUE: list of all genes
sDEG_dic = {}       # KEY: disease name  VALUE: list of sDEGs
df_dic = {}         # KEYS: disease      VALUE: its DE table

dis_list = list(disease_dirs) # disease list
for k in range(0,len(disease_dirs)):
    dis = dis_list[k]
    df = pd.read_table(results_dir+dis+'/'+dis+'_DEGs.txt') 
    filtered = df[df['adj.P.Val'] <= 0.05]
    allgenes_dic[dis] = set(df['symbol'])
    sDEG_dic[dis] = set(filtered['symbol'])
    df_dic[dis] = df
    

# Adding to the dictionaries the metapatient's information information --------------------------
#dis = 'AdenomatousPolyps'
m_dis_list = list(m_disease_dirs) # disease list
for k in range(0,len(m_disease_dirs)):
    dis = dis_list[k]
    try:
        m_dis_files = os.listdir(m_results_dir+dis)
    except:
        print("\nWARNING: There are no metapatients for "+dis)
        continue
    m_metapatients = [i for i in m_dis_files if '_DEGs_' in i] # Selecting files that correspond to sDEGs metapatient tables
    m_metapatients = [i for i in m_dis_files if ('_' + str(crandomization)+ '_') in i] # Select the meta-patients that correspond to the current randomization
    m_metapatients.sort()
    for k2 in range(0,len(m_metapatients)):
        df = pd.read_table(m_results_dir+dis+'/'+m_metapatients[k2]) 
        filtered = df[df['adj.P.Val'] <= 0.05]
        metap_name = dis+"_"+str((k2+1))
        allgenes_dic[metap_name] = set(df['symbol'])
        sDEG_dic[metap_name] = set(filtered['symbol'])
        #dm_dic[dis] = dm_table[['genes',dis]]     # UNCOMMENT FOR DM ANALYSIS 
        df_dic[metap_name] = df
        dis_list.append(metap_name)
        
# Updating values
n_diseases = len(dis_list)   
print(dis_list)

# Obtain the UNION of genes and sDEGs for ALL diseases
allgene_union = set(chain(*allgenes_dic.values())) # 20020
sgene_union = set(chain(*sDEG_dic.values())) # 17063

print(df_dic)

# Distance functions ---------------------------------------------------------
valid_distance_funcs = ['comparable_spearman_distance',
                        'pairwise_spearman_distance',
                        'pairwise_union_spearman_distance']


    
# Build df will all the sDEGs / all genes with a 0 logFC - to complete expression values of a given disease
zero_sgenes_df = pd.DataFrame(dict(zip(list(sgene_union),list(np.zeros(len(sgene_union))))),index=[0])
zero_sgenes_df = zero_sgenes_df.transpose() ;  zero_sgenes_df.reset_index(inplace=True) ; 
zero_sgenes_df.rename(columns={'index':'symbol',0:'logFC'}, inplace = True)
zero_allgenes_df = pd.DataFrame(dict(zip(list(allgene_union),list(np.zeros(len(allgene_union))))),index=[0])
zero_allgenes_df = zero_allgenes_df.transpose() ;  zero_allgenes_df.reset_index(inplace=True) ; 
zero_allgenes_df.rename(columns={'index':'symbol',0:'logFC'}, inplace = True)

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
       
#       print(spearmanr(list(expr1['logFC']), list(expr2['logFC'])))
       
       if(len(expr1['symbol']) > 1):
           if((np.var(list(expr1['logFC'])) != 0) and (np.var(list(expr2['logFC'])) != 0)):
               val,pv = spearmanr(list(expr1['logFC']), list(expr2['logFC']))
           else:
               val,pv = (np.nan,np.nan)
       else:
           val,pv = (np.nan,np.nan)
           
       
       return(val,pv)
       
# Create the distance matrix ------------------------------------------------- 
def build_distance_matrix(dis_list,distance_func, genes='sDEGs', sort=True, save_results=True):
    ''''''
    # Create an empty matrix
#    matrix = np.zeros((n_diseases,n_diseases))
    
    if ((distance_func.__name__ == 'pearson_distance') or (distance_func.__name__ == 'spearman_distance') or (distance_func.__name__ == 'comparable_spearman_distance') or (distance_func.__name__ == 'pairwise_union_spearman_distance') or (distance_func.__name__ == 'pairwise_spearman_distance')):
        outdf = pd.DataFrame(columns=['Dis1','Dis2','Distance','pvalue']) #output df
#        pv_matrix = np.zeros((n_diseases,n_diseases))
       
        for k1 in range(0,len(dis_list)):
            for k2 in range(k1+1,len(dis_list)):
                try:
                    print(dis_list[k1] + '  ' + dis_list[k2])
                    cdis,cpv = distance_func(dis_list[k1],dis_list[k2],genes=genes)
#                    print(dis_list[k1]+" "+dis_list[k2]+" --> "+str('%.4f' %cdis)+"  pv= "+('%.4f' %cpv))
#                    matrix[k1, k2] = cdis
#                    pv_matrix[k1, k2] = cpv
                    outdf = outdf.append({'Dis1':dis_list[k1],'Dis2':dis_list[k2],'Distance':cdis,'pvalue':cpv},ignore_index=True)
                except:
                    print(dis_list[k1] + ' ' + dis_list[k2] + ' error')
                
        if sort:
            outdf = outdf.sort_values(by=['Distance','pvalue'], ascending=False, inplace=False)
        if save_results:
            outpath = outpath = 'distances/'+ distance_func.__name__ + "_" + genes + "_" + str(crandomization) + '.csv'
            outdf.to_csv(outpath, index=False)
        return outdf
        
        
    elif (distance_func.__name__ in valid_distance_funcs):
        outdf = pd.DataFrame(columns=['Dis1','Dis2','Distance']) #output df
        for k1 in range(0,len(dis_list)):
            for k2 in range(k1+1,len(dis_list)):
                cdis = distance_func(dis_list[k1],dis_list[k2],genes=genes)
#                print(dis_list[k1]+" "+dis_list[k2]+" --> "+str(cdis))
#                matrix[k1, k2] = cdis
                outdf = outdf.append({'Dis1':dis_list[k1],'Dis2':dis_list[k2],'Distance':cdis},ignore_index=True)
                
        if sort:
            outdf = outdf.sort_values(by='Distance', ascending=False, inplace=False, na_position='last')
        if save_results:
            outpath = 'distances/'+ distance_func.__name__ + "_" + genes + "_" + str(crandomization) + '.csv'
            outdf.to_csv(outpath, index=False)
        return outdf
    
    else:
        raise ValueError("Non valid distance_func")
        
outdf = build_distance_matrix(dis_list, pairwise_union_spearman_distance) 

# Obtain the network from the distances by correcting for multiple testing (FDR) --------------------------------------

# Dealing with NaNs
pvals = outdf['pvalue']
#pvals[1:30]
mask = np.isfinite(pvals)    # True if a number, False if Infinite or NaN
pval_corrected = np.empty(pvals.shape)
pval_corrected.fill(np.nan) 
pval_corrected[mask] = fdrcorrection(pvals[mask],0.05, is_sorted=False)[1]
outdf['adj.pvalue'] = pval_corrected

final_network = outdf[outdf['adj.pvalue'] < 0.05]
final_network.shape
outpath = 'defined_networks/'+ 'network_pairwise_union_spearman_distance' + "_" + 'sDEGs' + "_" + str(crandomization) + '.csv'
final_network.to_csv(outpath, index=False)
    
