#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 18:40:20 2022

@author: beatriz
"""

# Network backbone

import networkx as nx
import distanceclosure as dc
import pandas as pd
import matplotlib.pyplot as plt

def prepare_network_for_backbone(net, epsilon):
    #net = dsn
    pos = net[net.Distance > 0]; print(pos.shape)
    neg = net[net.Distance < 0]; print(neg.shape)
    
    dfl = (pos, neg)
    odf = []
    for df in dfl:
        df = pd.DataFrame({'source':df.Dis1,
                             'target':df.Dis2,
                             'distance':(1 - abs(df.Distance) + epsilon)   
        })
        odf.append(df)
        
    return(odf)
    

def get_network_backbone(df, metric='metric'):
    #df = odf[1]
    G = nx.from_pandas_edgelist(df, edge_attr=True)
    G = G.to_undirected()
    print("Network: " + str(G.number_of_edges()))
    backbone = dc.distance_closure(G, kind=metric, weight='distance', only_backbone=True)
#    print("Backbone: " + str(backbone.number_of_edges()))
    nx.draw(backbone, with_labels=True)
    backbonedf = nx.to_pandas_edgelist(backbone)
    if metric == 'metric':
        filtered = backbonedf[backbonedf.is_metric == True]
    else:
        filtered = backbonedf[backbonedf.is_ultrametric == True]
    nbackbone = filtered.shape[0]
    print("Backbone: " + str(nbackbone))
    print("Backbone: " + str(100*(nbackbone/G.number_of_edges())))
    return(backbonedf, nbackbone)
    
def merge_and_save_backbones(dfl, filename):
#    dfl = dsn_backbones
#    filename = "shiny_network_app/additional_data/backbones/dsn_metric.txt"
    
    df = pd.concat([dfl[0], dfl[1]])
    df.to_csv(filename, sep="\t", index=False)

dsn = pd.read_csv("shiny_network_app/final_pairwise_union_spearman_distance_sDEGs_network.txt", sep="\t")
ssn = pd.read_csv("shiny_network_app/final_metap_dis_pairwise_union_spearman_distance_sDEGs_network.txt", sep="\t")

dim_dsn = dsn.shape[0]
dim_ssn = ssn.shape[0]

epsilon= 0.000001

# DSN
print("\nDSN:")
dsnl = prepare_network_for_backbone(dsn, epsilon)

# metric
print("\nMetric:")
dsn_backbones = []
nbackbone = 0
for dsn in dsnl:
    backbone, cnbackbone = get_network_backbone(dsn)
    nbackbone += cnbackbone
    dsn_backbones.append(backbone)

filename = "shiny_network_app/additional_data/backbones/dsn_metric.txt" 
merge_and_save_backbones(dsn_backbones, filename)
global_backbone_perc = (nbackbone / dim_dsn) * 100; print("\nGlobal backbone %: " + str(global_backbone_perc))

# ultrametric
print("\nUltrametric:")
dsn_backbones2 = []
nbackbone = 0
for dsn in dsnl:
    backbone, cnbackbone  = get_network_backbone(dsn, metric='ultrametric')
    nbackbone += cnbackbone
    dsn_backbones2.append(backbone)
    
filename = "shiny_network_app/additional_data/backbones/dsn_ultrametric.txt" 
merge_and_save_backbones(dsn_backbones2, filename)
global_backbone_perc = (nbackbone / dim_dsn) * 100; print("\nGlobal backbone %: " + str(global_backbone_perc))
    
    
# SSN
print("\nSSN:")
ssnl = prepare_network_for_backbone(ssn, epsilon)

# metric
print("\nMetric:")
ssn_backbones = []
nbackbone = 0
for ssn in ssnl:
    backbone, cnbackbone  = get_network_backbone(ssn)
    nbackbone += cnbackbone
    ssn_backbones.append(backbone)

filename = "shiny_network_app/additional_data/backbones/ssn_metric.txt" 
merge_and_save_backbones(ssn_backbones, filename)
global_backbone_perc = (nbackbone / dim_ssn) * 100; print("\nGlobal backbone %: " + str(global_backbone_perc))

# ultrametric
print("\nUltrametric:")
ssn_backbones2 = []
nbackbone = 0
for ssn in ssnl:
    backbone, cnbackbone  = get_network_backbone(ssn, metric='ultrametric')
    nbackbone += cnbackbone
    ssn_backbones2.append(backbone)

filename = "shiny_network_app/additional_data/backbones/ssn_ultrametric.txt" 
merge_and_save_backbones(ssn_backbones2, filename)
global_backbone_perc = (nbackbone / dim_ssn) * 100; print("\nGlobal backbone %: " + str(global_backbone_perc))
    
