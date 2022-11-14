#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 10:48:36 2022

This script is used to analysze and plot subsets of the AIBS RNA-seq database (https://portal.brain-map.org/atlases-and-data/rnaseq).
MTG - SMART-SEQ (2019) seperate matrices (exons and introns) is used to sub-select expression values
for specific genes of interest per transcriptomic type within MTG. 

@author: annagalakhova & stand
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#%% LOAD AND PROCESS METADATA FILES 
#load the metadata to select samples 
samples_metadata = pd.read_csv(r"path\human_MTG_2018-06-14_samples-columns.csv")
#select samples based on layers  
layers = ['L2', 'L3']
samples_metadata_L23 = samples_metadata[samples_metadata['brain_subregion'].isin(layers)] 

#%%SELECT TRANSCRIPTOMIC TYPES
t_types = ['Exc L2-3 LINC00507 FREM3','Exc L2 LAMP5 LTK','Exc L2-4 LINC00507 GLP2R','Exc L3-4 RORB CARM1P1','Exc L3-5 RORB COL22A1']
samples_metadata_L23_types =samples_metadata_L23[samples_metadata_L23['cluster'].isin(t_types)]
#change the class label to Class, since class is recoginzed as a python command
samples_metadata_L23_types['class'] = samples_metadata_L23_types.rename(columns = {'class':'Class'}, inplace = True)

#%% LOAD GENES METADATA AND THE GENE SET OF INTEREST 
#load the gene metadata, containing information about all the genes. 
gene_metadata = pd.read_csv(r"path\human_MTG_2018-06-14_genes-rows.csv")

#load in the geneset from text_file NOTE: add quotechar if apostrophe's are present in the original text file (i.e. 'gene' instead of gene)
geneset = pd.read_csv(r"path\geneset.txt", names=['gene'], header=None , quotechar="'")
geneset = geneset.gene.drop_duplicates()
geneset = pd.DataFrame(geneset)
#creating a variable that contains all the gene names
all_genes = gene_metadata.gene

#%% LOAD IN THE DATA
#get expression data of genes of interest from exons and introns
skiprows_genes = all_genes.loc[~all_genes.isin(geneset.gene)].index+1
#samples to load from the matrix 
usecols = list(samples_metadata_L23_types.sample_name)

#%%run this to check what genes are not present in the original matrices
keeprows =  all_genes.loc[all_genes.isin(geneset.gene)].index+1
genes_not_pres = pd.DataFrame(all_genes.loc[all_genes.isin(geneset.gene)])
excluded = geneset.gene.loc[~geneset.gene.isin(genes_not_pres.gene)]

#%%LOAD A FILTERED MATRIX BASED ON GENE SET OF INTEREST
#load in the data from exons
exons = pd.read_csv(r"path/exon_mtg_log10_CPM.csv", 
                        skiprows=skiprows_genes, usecols=usecols)
exons.loc[:,'gene']=geneset
exons=exons.set_index('gene').T
#load in the data from introns
introns = pd.read_csv(r"path/intron_mtg_log10_CPM.csv",
                        skiprows=skiprows_genes, usecols=usecols)
introns.loc[:,'gene']=geneset
introns=introns.set_index('gene').T
#add exons + introns
data_log10CPM = exons + introns


#%%CLEAN THE DATA FROM NON-EXPRESSING GENES
#crucial in getting max per column values, transposing columns containg only zero will fill in NaN's, replace back to zeros 
data_log10CPM = data_log10CPM.replace(np.nan,0)
#obtaining the maximum expression of one gene across all the cells
max_per_column = data_log10CPM.loc[:,:].max(axis=0)
#loop over the maximum values per gene to identify zero expressing genes, and append these genes in a list. 
excluded_genes = []
for i in range(0,len(max_per_column)):
    temp = max_per_column[i]
    if temp == 0:
        excluded_genes.append(max_per_column.index[i])
    del temp 
#delete genes that are cut off (above) from the dataset   
data_log10CPM = data_log10CPM[data_log10CPM.columns.drop(excluded_genes)]

#%% SAVE OR LOAD PRE-PROCESSED DATASETS
data_log10CPM.to_csv(r"path/example_name.csv")


#%%ADDING AVERAGE EXPRESSION VALUES PER CELL
data_log10CPM['avg'] = data_log10CPM.mean(axis=1)

#%% IN CASE NEEDED TO NORMALIZE TO HOUSEKEEPING GENES, LOAD IN THE HOUSEKEEPING GENES DATASET
# load the HK values to calculate ratio between HK and gene set of interest  
HK_log10CPM = pd.read_csv(r"path/example_name.csv")
#get average HK column 
average_HK = HK_log10CPM.avg

data_log10CPM['avg_correct'] =  list(average_HK)/ data_log10CPM['avg'] 

#%% ADDING METADATA TO THE DATASET
#add sample name, region and class information to the dataframe 
sample_names = list(data_log10CPM.index)
samples_metadata_indexed = samples_metadata_L23_types.loc[samples_metadata_L23_types['sample_name'].isin(sample_names)]
#add layer, t-type and class info to the dataframe 
data_log10CPM['brain_subregion'] = list(samples_metadata_indexed.brain_subregion)
data_log10CPM['cluster'] = list(samples_metadata_indexed.cluster)
data_log10CPM['Class'] = list(samples_metadata_indexed.Class)


#%% SPLIT THE DATA BASED ON T-TYPE
LTK_avg = data_log10CPM.loc[data_log10CPM.cluster.values == 'Exc L2 LAMP5 LTK']
GLP2R_avg =  data_log10CPM.loc[data_log10CPM.cluster.values == 'Exc L2-4 LINC00507 GLP2R']
FREM3_sup_avg =  data_log10CPM.loc[(data_log10CPM.cluster.values == 'Exc L2-3 LINC00507 FREM3') & (data_log10CPM.brain_subregion.values == 'L2')]
FREM3_deep_avg = data_log10CPM.loc[(data_log10CPM.cluster.values == 'Exc L2-3 LINC00507 FREM3') & (data_log10CPM.brain_subregion.values == 'L3')]
CARM1P1_avg = data_log10CPM.loc[data_log10CPM.cluster.values =='Exc L3-4 RORB CARM1P1']
COL22A1_avg = data_log10CPM.loc[data_log10CPM.cluster.values =='Exc L3-5 RORB COL22A1']

#%%plot the data for t_types

ax = sns.violinplot(data=[LTK_avg.avg_correct,  FREM3_sup_avg.avg_correct, FREM3_deep_avg.avg_correct, GLP2R_avg.avg_correct, CARM1P1_avg.avg_correct, COL22A1_avg.avg_correct], inner = 'quartile',scale="width", linewidth=4,
                   alpha=0.8, figsize=(10,16))
ax.set_xticklabels(['LTK', "L2-FREM3", "L3-FREM3", "GLP2R","CARM1P1", "COL22A1"])
ax.set_xlabel('cell type')
ax.set_ylabel("Ratio HK/HAR Mean Log10-(CPM+1)")
#ax.set_ylim(0.7,1.5)
for l in ax.lines:
    l.set_linestyle('--')
    l.set_linewidth(0)
    l.set_color('red')
    l.set_alpha(0.8)
for l in ax.lines[1::3]:
    l.set_linestyle('-')
    l.set_linewidth(3)
    l.set_color('black')
    l.set_alpha(0.8)
    
plt.savefig("ttypes_HAR_HK_ratio.eps", format='eps')
plt.savefig("ttypes_HAR_HK_ratio.png", format ='png')






















