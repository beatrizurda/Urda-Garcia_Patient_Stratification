#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 15:55:08 2020

@author: beatriz
"""

import pandas as pd
import glob
# import re

# To test
#cmetapdir = 'Metapatients/with_entire_count_matrix/FE_results/BreastCancer_DEGs_allgenes3_1.GseaPreranked.1589148231433'


class Deregulated_pathways:
    def __init__(self,df_up, df_down):
        self.up = df_up[df_up['adj.P.Val'] <= 0.05]
        self.down = df_down[df_down['adj.P.Val'] <= 0.05]
    

class Metapatient:
    def __init__(self, name):
        self.name = name
        self.disease = name.split("_")[0]
        self.number = name.split("_")[1]
        self.sDEGs = []
        self.sDVGs = []
        self.DE_pathways = []
        self.DV_pathways = []
        self.sDEGs_up = []
        self.sDEGs_down = []
        self.DE_pathways_up = []
        self.DE_pathways_down = []
        
    def get_sDEGs(self):
        if not self.sDEGs:   # if it is empty
            cdir = "Metapatients/with_entire_count_matrix/DEAnalysis/"+self.disease+"/"
            df_degs = pd.read_table(glob.glob(cdir+self.disease+"_DEGs_"+"*_"+str(self.number)+".txt")[0])
            #print(df_degs.dtypes)
            df_degs = df_degs[df_degs['adj.P.Val'] < 0.05]
            df_degs = df_degs.sort_values("logFC", ascending=False)  # sort it by logFC
            self.sDEGs = df_degs
            self.sDEGs_up = df_degs[df_degs['logFC'] > 0]
            self.sDEGs_down = df_degs[df_degs['logFC'] < 0].sort_values("logFC")
            return(self.sDEGs_up, self.sDEGs_down)
        else:
            return(self.sDEGs_up, self.sDEGs_down)
    
    def sDVGs(self):
        pass
    
    def get_DE_pathways(self):
         if not self.DE_pathways:   # if it is empty
            cdir = "Metapatients/with_entire_count_matrix/FE_results/"
            # list.files(path=input_dir, pattern = paste(dis,"_",sep=""))
            cmetapdir = glob.glob(cdir+self.disease+"_*_"+str(self.number)+"*")[0]; print(cmetapdir)
            # return(cmetapdir)
            df_up = pd.read_table(glob.glob(cmetapdir+"/"+"gsea_report_for_na_pos_*.xls")[0])
            df_up.columns = df_up.columns.str.replace(' ', '.')
            #df_up['FDR.q-val']
            df_up = df_up[df_up['FDR.q-val'] <= 0.05]
            # file_up = pd.read_table(glob.glob('Metapatients/with_entire_count_matrix/FE_results/BreastCancer_DEGs_allgenes3_1.GseaPreranked.1589148231433'+"/"+"gsea_report_for_na_pos_*.xls")) 
            df_down = pd.read_table(glob.glob(cmetapdir+"/"+"gsea_report_for_na_neg_*.xls")[0])
            df_down.columns = df_down.columns.str.replace(' ', '.')
            df_down = df_down[df_down['FDR.q-val'] <= 0.05]
            self.DE_pathways_up = set(df_up['NAME'])
            self.DE_pathways_down = set(df_down['NAME'])
            self.DE_pathways = (df_up, df_down)
            
            return(df_up, df_down)
         else:
            return self.DE_pathways
    
    def DV_pathways(self):
        pass
    
class Disease:
    def __init__(self, name):
        self.name = name
        self.n_metapatients = []
        self.sDEGs = []
        self.sDVGs = []
        self.DE_pathways = []
        self.DV_pathways = []
        self.sDEGs_up = []
        self.sDEGs_down = []
        self.DE_pathways_up = []
        self.DE_pathways_down = []
        
    def get_sDEGs(self):
        if not self.sDEGs:   # if it is empty
            cdir = "DL_RESULTS_GREIN/"+self.name+"/"
            df_degs = pd.read_table(glob.glob(cdir+"/"+self.name+"_DEGs.txt")[0])
            print(df_degs.dtypes)
            df_degs = df_degs[df_degs['adj.P.Val'] < 0.05]
            
            df_degs = df_degs.sort_values("logFC", ascending=False)  # sort it by logFC
            print(df_degs)
            self.sDEGs = df_degs
            self.sDEGs_up = df_degs[df_degs['logFC'] > 0]
            self.sDEGs_down = df_degs[df_degs['logFC'] < 0].sort_values("logFC")
            return(self.sDEGs_up, self.sDEGs_down)
        else:
            return(self.sDEGs_up, self.sDEGs_down)
    
    def sDVGs(self):
        pass
    
    def get_DE_pathways(self):
        if not self.DE_pathways:   # if it is empty
            cdir = "FE_results/Pathways_DEGs_diseases/"
            cdisdir = glob.glob(cdir+self.name+"_*")[0]; print(cdisdir)
            
            df_up = pd.read_table(glob.glob(cdisdir+"/"+"gsea_report_for_na_pos_*.xls")[0])
            df_up.columns = df_up.columns.str.replace(' ', '.')
            #df_up['FDR.q-val']
            df_up = df_up[df_up['FDR.q-val'] <= 0.05]
            # file_up = pd.read_table(glob.glob('Metapatients/with_entire_count_matrix/FE_results/BreastCancer_DEGs_allgenes3_1.GseaPreranked.1589148231433'+"/"+"gsea_report_for_na_pos_*.xls")) 
            df_down = pd.read_table(glob.glob(cdisdir +"/"+"gsea_report_for_na_neg_*.xls")[0])
            df_down.columns = df_down.columns.str.replace(' ', '.')
            df_down = df_down[df_down['FDR.q-val'] <= 0.05]
            self.DE_pathways_up = set(df_up['NAME'])
            self.DE_pathways_down = set(df_down['NAME'])
            self.DE_pathways = (df_up, df_down)
            return(df_up, df_down)
        else:
            return self.DE_pathways
            
    def DV_pathways(self):
        pass
        
    
    
####### OVERALL DYSREGULATED PATHWAYS ###########
reactome_parents = pd.read_table('Reactome/Reactome_parents.txt') 
metadata = pd.read_table('new_disease_metadata_final_names.txt') 
icd9df = metadata[['disease_name','icd9','disease_cat']]
dis_list = metadata['disease_name']  
ndis = len(dis_list) 

# First create a disctionary with the DISEASE OBJECTS:
# Key: disease_name | Values: disease object (quick access to pathwats_up, pathwats_down)
dis_dic = {}
for disname in dis_list:
    cdis = Disease(disname)
    cdis.get_DE_pathways()
    dis_dic[disname] = cdis
    
# Ejermplo de uso:
#dis_dic['Glioblastoma'].DE_pathways_up
#dis_dic['BorreliaBurdorferiInfection'].DE_pathways_down
    

common_pathways = []
for k1 in range(ndis):
    dis1 = dis_list[k1]
    for k2 in range(k1+1,ndis):
        dis2 = dis_list[k2]
        cup = dis_dic[dis1].DE_pathways_up.union(dis_dic[dis2].DE_pathways_up)
        cdown = dis_dic[dis1].DE_pathways_down.union(dis_dic[dis2].DE_pathways_down)
        if((len(cup) == 0) or (len(cdown) == 0)):
            print([dis1,dis2,cup, cdown])
        common_pathways.append([dis1,dis2,cup, cdown])
        
common_pathwaysdf = pd.DataFrame(common_pathways, columns=['dis1','dis2','common_up', 'common_down'])

common_pathwaysdf.to_csv("common_pathways_disease_level.csv", index=False)

# Add the disease category of Dis1 and Dis1 
pass


### MERGE with ICD9 codes and disease categories
common_pathwaysdf = pd.merge(common_pathwaysdf, icd9df, how='left', left_on='dis1', right_on='disease_name')
common_pathwaysdf = common_pathwaysdf.drop(['disease_name'], axis=1)
common_pathwaysdf.columns = ['dis1', 'dis2', 'common_up', 'common_down','dis1_icd9','dis1_category']

common_pathwaysdf = pd.merge(common_pathwaysdf, icd9df, how='left', left_on='dis2', right_on='disease_name')
common_pathwaysdf = common_pathwaysdf.drop(['disease_name'], axis=1)
common_pathwaysdf.columns = ['dis1', 'dis2', 'common_up', 'common_down','dis1_icd9','dis1_category', 'dis2_icd9','dis2_category']
common_pathwaysdf.dtypes

### SORT INTERACTIONS smaller to bigger
#for index,row in common_pathwaysdf.head(n=2).iterrows():
#    print(index,row)
    
for index,row in common_pathwaysdf.iterrows():
#    print(index,row)
#    print(row['dis1_icd9'], row['dis2_icd9'])
    if(row['dis1_icd9'] > row['dis2_icd9']):
        common_pathwaysdf.loc[index, 'dis1_icd9'], common_pathwaysdf.loc[index, 'dis2_icd9'] = row['dis2_icd9'], row['dis1_icd9']
        common_pathwaysdf.loc[index, 'dis1'], common_pathwaysdf.loc[index, 'dis2'] = row['dis2'], row['dis1']
        common_pathwaysdf.loc[index, 'dis1_category'], common_pathwaysdf.loc[index, 'dis2_category'] = row['dis2_category'], row['dis1_category']
    
        
# Checking it's good.     
for index,row in common_pathwaysdf.head(n=5).iterrows():
#    print(index,row)
    print(row['dis1_icd9'], row['dis2_icd9'])
    
### FILTERING INTERACITIONS IN EPIDEMIOLOGY AND NOT IN EPIDEMIOLOGY
common_pathwaysdf['interactions'] =  common_pathwaysdf['dis1_icd9'].astype(str) + '_' + common_pathwaysdf['dis2_icd9'].astype(str)

# Divide in TWO: IN EPIDEMIOLOGY and NOT IN EPIDEMIOLOGY
dsn_path = 'Network_building/Defined_networks/pairwise_union_spearman_distance_sDEGs_network.txt'
dsn = pd.read_table(dsn_path) 
dsn_pos = dsn.loc[dsn['Distance'] > 0]

# PROCESSED DSN. At the disease level transformed into icd9 codes.
# Only positive interactions, in and not in epidemiology. 
processed_pos_dsn  = pd.read_table('Network_building/Overlapping_results/Shuffling_labels/pairwise_union_spearman_distance_sDEGs_pos_consistent_B_TRUE_final_network.txt')
processed_pos_dsn['interactions'] = processed_pos_dsn['Dis1'].astype(str) + '_' + processed_pos_dsn['Dis2'].astype(str)


# PROCESSED BARABASI NETWORK. icd9 codes, sorted, in DSN.
processed_barabasi = pd.read_table('Network_building/Overlapping_results/Shuffling_labels/pairwise_union_spearman_distance_sDEGs_pos_consistent_B_TRUE_final_barabasi.txt')
processed_barabasi['interactions'] = processed_barabasi['Dis1'].astype(str) + '_' + processed_barabasi['Dis2'].astype(str)

processed_pos_dsn.shape # 304 = 130 + 174
in_epidem = processed_pos_dsn.loc[processed_pos_dsn['interactions'].isin(processed_barabasi['interactions'])]
not_epidem = processed_pos_dsn.loc[~processed_pos_dsn['interactions'].isin(processed_barabasi['interactions'])]

in_epidem.shape # 130     ---> 129
not_epidem.shape # 174    ---> 175

mean_corr1 = sum(in_epidem['Distance'])/len(in_epidem['Distance'])
mean_corr2 = sum(not_epidem['Distance'])/len(not_epidem['Distance'])
mean_corr1/mean_corr2

from scipy import stats
tdsn, pdsn = stats.ttest_ind(in_epidem['Distance'], not_epidem['Distance'])


#### ADD THE PATHWAY INFORMATION
# That is, from common_pathwaysdf filter the rows that are in_epidem.
pathways_in_epidem = common_pathwaysdf.loc[common_pathwaysdf['interactions'].isin(in_epidem['interactions'])]
pathways_not_epidem = common_pathwaysdf.loc[common_pathwaysdf['interactions'].isin(not_epidem['interactions'])]


pathways_in_epidem.shape # 146       ---> 145   ICD: 44, 88--> 154
pathways_not_epidem.shape # 184      ---> 185   ICD: 44, 88--> 197

pathways_in_epidem['n_common_up'] = list(map(len, pathways_in_epidem['common_up']))
pathways_in_epidem['n_common_down'] = list(map(len, pathways_in_epidem['common_down']))

pathways_not_epidem['n_common_up'] = list(map(len, pathways_not_epidem['common_up']))
pathways_not_epidem['n_common_down'] = list(map(len, pathways_not_epidem['common_down']))

pathways_in_epidem.to_csv("PLOTS/Pathway_counts/pathways_dsn_in_epidemiology.csv",  sep="\t", index=False)
pathways_not_epidem.to_csv("PLOTS/Pathway_counts/pathways_dsn_not_in_epidemiology.csv",  sep="\t", index=False)


mean_up_epidem = sum(pathways_in_epidem['n_common_up'])/(pathways_in_epidem.shape[0]); mean_up_epidem
mean_down_epidem = sum(pathways_in_epidem['n_common_down'])/(pathways_in_epidem.shape[0]); mean_down_epidem

mean_up_not_epidem = sum(pathways_not_epidem['n_common_up'])/(pathways_not_epidem.shape[0]); mean_up_not_epidem
mean_down_not_epidem = sum(pathways_not_epidem['n_common_down'])/(pathways_not_epidem.shape[0]); mean_down_not_epidem


#### COUNT the number of interactions between each disease category pair
#pathways_in_epidem['discat_pair'] = pathways_in_epidem['dis1_category'] + '_' + pathways_in_epidem['dis2_category']
#pathways_not_epidem['discat_pair'] = pathways_not_epidem['dis1_category'] + '_' + pathways_not_epidem['dis2_category'] 


# COUNT THE NUMBER OF OCCURRENCES FOR EACH PATHWAY

# Testing the 2 lines below: creating a set from a list of sets
#hey = [{1,2}, {2,3,4,8}, {'a','b',8,9}]
#b=set().union(*hey)
#b

all_up = set().union(*common_pathwaysdf['common_up']); len(all_up) # 3326
all_down = set().union(*common_pathwaysdf['common_down']); len(all_down) # 3077

all_pathways = all_up.union(all_down); len(all_pathways) # 4110

def get_pathways_count(df=pathways_in_epidem):
    '''
    Returns a 
    df = pathways_in_epidem
    '''
    outputdf = []
    for pathway in all_pathways:
        up = []
        down = []
        up_count = 0; down_count = 0;
        
        for index,row in df.iterrows():  # For each interaction
            if(pathway in row['common_up']):
                up.append(row['dis1_category']+'_'+row['dis2_category'])
                up_count += 1
            if(pathway in row['common_down']):
                down.append(row['dis1_category']+'_'+row['dis2_category'])
                down_count += 1
                
#        print((pathway, up_count, down_count, up, down))
        outputdf.append([pathway, up_count, down_count, up, down])
        
    outputdf = pd.DataFrame(outputdf, columns=['pathways', 'up_count', 'down_count', 'up', 'down'])    
    return outputdf


output_epidem = get_pathways_count(pathways_in_epidem)
output_not_epidem = get_pathways_count(pathways_not_epidem)

### ADD PATHWAY PARENT CATEGORY
output_epidem = pd.merge(output_epidem, reactome_parents, how='left', left_on='pathways', right_on='GSEA_pathway')
output_epidem = output_epidem.drop(['GSEA_pathway'], axis=1)

output_not_epidem = pd.merge(output_not_epidem, reactome_parents, how='left', left_on='pathways', right_on='GSEA_pathway')
output_not_epidem = output_not_epidem.drop(['GSEA_pathway'], axis=1)


output_epidem.to_csv("pathways_count_dsn_in_epidemiology.csv", sep="\t", index=False)
output_not_epidem.to_csv("pathways_count_dsn_not_in_epidemiology.csv",  sep="\t", index=False)
                
        
epidem_reactome = output_epidem.loc[output_epidem['pathways'].isin(reactome_parents['GSEA_pathway'])]
not_epidem_reactome = output_not_epidem.loc[output_not_epidem['pathways'].isin(reactome_parents['GSEA_pathway'])]
# Shape: 801 Reactome pathways. 

         
parents = list(set(reactome_parents['Parent']))

### BEFORE - without the interactions that are involved
#parent_level_count = []
#for parent in parents:
##    parent = parents[0]   # Comment
#    crow = [parent]                                      
#    for cdf in [epidem_reactome, not_epidem_reactome]:
##        cdf = epidem_reactome                                  # Comment
#        cdf = cdf.loc[cdf['Parent'] == parent]
#        cup = sum(cdf['up_count']); cdown = sum(cdf['down_count'])
#        crow.extend([cup, cdown])
#    parent_level_count.append(crow)
#        
#parent_level_count = pd.DataFrame(parent_level_count, columns=['Parent','epidem_up','epidem_down', 'not_epidem_up','not_epidem_down']) 
#parent_level_count.to_csv("PLOTS/Pathway_counts/parent_level_counts.csv",  sep="\t", index=False)  

from itertools import chain

parent_level_count = []
for parent in parents:
#    parent = parents[0]   # Comment
    crow = [parent]                                      
    for cdf in [epidem_reactome, not_epidem_reactome]:
#        cdf = epidem_reactome                                  # Comment
        
        cdf = cdf.loc[cdf['Parent'] == parent]
        up_count = sum(cdf['up_count']); down_count = sum(cdf['down_count'])
        up=list(chain(*cdf['up']))   # add the elements of an indefinite number of lists
        down=list(chain(*cdf['down']))
        crow.extend([up_count, down_count, up, down])
    parent_level_count.append(crow)
    
parent_level_count = pd.DataFrame(parent_level_count, columns=['Parent',
                                                               'epidem_up','epidem_down', 
                                                               'epidem_up_discat','epidem_down_discat',
                                                               'not_epidem_up','not_epidem_down',
                                                               'not_epidem_up_discat','not_epidem_down_discat']) 
    
parent_level_count.to_csv("PLOTS/Pathway_counts/parent_level_counts.csv",  sep="\t", index=False)  


interactions_by_parents = []
for parent in parents:
#    parent = parents[0]   # Comment
    crow = [parent]
    cpathways = list(reactome_parents.loc[reactome_parents['Parent'] == parent]['GSEA_pathway'])
    for cdf in [pathways_in_epidem, pathways_not_epidem]:
#        cdf = pathways_in_epidem   # Comment
        
        up_count = 0; down_count = 0
        up = []; down = []
        for index, row in cdf.iterrows(): # for each line in the df
            if(len(set(cpathways).intersection(row['common_up'])) > 0):
                up_count += 1
                up.append(row['dis1_category']+'_'+row['dis2_category'])
            if(len(set(cpathways).intersection(row['common_down'])) > 0):
                down_count += 1
                down.append(row['dis1_category']+'_'+row['dis2_category'])
                
        
        crow.extend([up_count, down_count, up, down])
    
    
    interactions_by_parents.append(crow)            

interactions_by_parentsdf = pd.DataFrame(interactions_by_parents, columns=['Parent', 
                                                                           'epidem_up','epidem_down', 
                                                                           'epidem_up_discat','epidem_down_discat',
                                                                           'not_epidem_up','not_epidem_down',
                                                                           'not_epidem_up_discat','not_epidem_down_discat'])       
    

interactions_by_parentsdf.to_csv("PLOTS/Pathway_counts/parent_level_counts_explained_interactions.csv",  sep="\t", index=False)


    
    
    
###### HUNGTINTON DISEASE WITH LIVER, LUNG AND BREAST CANCER ###########
# HTT ENSG00000197386
HD = Disease("HuntingtonsDisease"); HD_pathways = HD.get_DE_pathways(); HD_sdegsup, HD_sdegsdown = HD.get_sDEGs()
Liverc = Disease("LiverCancer"); Liverc_pathways = Liverc.get_DE_pathways(); Liverc_sdegsup, Liverc_sdegsdown  = Liverc.get_sDEGs()
Lungc = Disease("LungCancer"); Lungc_pathways = Lungc.get_DE_pathways(); Lungc_sdegsup, Lungc_sdegsdown = Lungc.get_sDEGs()   
Cll = Disease("ChronicLymphocyticLeukemia"); Cll_pathways = Cll.get_DE_pathways(); Cll_sdegsup,Cll_sdegsdown  = Cll.get_sDEGs()  
Breastcancer = Disease("BreastCancer"); Breastcancer_pathways = Breastcancer.get_DE_pathways(); Breastcancer_sdegsup, Breastcancer_sdegsdown = Breastcancer.get_sDEGs()  

"ENSG00000197386" in list(HD_sdegsup['symbol'])
"ENSG00000197386" in list(Liverc_sdegsup['symbol'])
"ENSG00000197386" in list(Lungc_sdegsup['symbol'])
"ENSG00000197386" in list(Cll_sdegsup['symbol'])
"ENSG00000197386" in list(Breastcancer_sdegsup['symbol'])

"ENSG00000197386" in list(HD_sdegsdown['symbol'])      # TRUE
"ENSG00000197386" in list(Liverc_sdegsdown['symbol'])  # F
"ENSG00000197386" in list(Lungc_sdegsdown['symbol'])   # TRUE
"ENSG00000197386" in list(Cll_sdegsdown['symbol'])     # F
"ENSG00000197386" in list(Breastcancer_sdegsdown['symbol']) # TRUE









