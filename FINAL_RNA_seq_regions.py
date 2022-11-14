#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 10:37:32 2022

This script is used to analyze and plot subsets of the AIBS RNA-seq database (https://portal.brain-map.org/atlases-and-data/rnaseq).
MULTIPLE CORTICAL AREAS - SMART-SEQ (2019) full matrix (exons+introns) is used to sub-select expression values
for specific genes of interest per region or cell classes. 

Script structure: 
19-63   selecting the data based on gene set of interest (If datasets needs to be re-used, data can be saved and loaded from lines 69 & 72 ) 
73-126  Further sub-selection, dividing the data into regions and cell-classes.
128-end Plotting the data 

@authors: annagalakhova & stand
"""


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#%% LOAD AND PROCESS METADATA FILES 
#load the metadata from the metadata file (cell-type-database cortical regions)
samples_metadata = pd.read_csv(r"path\metadata.csv")
#clear the metadata from samples with no class present.  
samples_metadata_cleared = samples_metadata[samples_metadata['class_label'].notnull()]
#%% LOAD THE GENE SET OF INTEREST 
#load in the geneset from text_file NOTE: add quotechar if apostrophe's are present in the original text file (i.e. 'gene' instead of gene)
geneset = pd.read_csv(r"path/geneset.txt", names=['gene'], header=None   , quotechar="'")
#incase there are duplicates present in the gene set, remove and update the gene set 
geneset = geneset.gene.drop_duplicates()
#convert geneset to a dataframe to improve compatibillity with pandas 
geneset = pd.DataFrame(geneset)

#%% LOAD A FILTERED MATRIX BASED ON GENE SET OF INTEREST 
# selecting expression data of the genes present in the gene set of interest (geneset)
skiprows_set1 = samples_metadata.loc[~samples_metadata.sample_name.isin(samples_metadata_cleared.sample_name)].index+1
#load in the data, NOTE. here we used pre-normalized log10 CPM data. 
data_log10CPM = pd.read_csv(r"path\matrix_log10_CPM.csv",
                  skiprows=skiprows_set1, usecols=geneset.gene)

#%%RUN THIS PART IF GENES (COLUMNS) ARE NOT FOUND IN THE MATRIX 
#running the cell above can result in an error that certain columns (i.e. genes) are not present in the matrix. 
#in that case copy and paste the gene names from the error under the variable deleted genes (see below)
deleted_genes =   ['LARGE1', 'GPAT3', 'SQOR', 'SVBP', 'FAM241A', 'CYTOR', 'MGC45800', 'LOC731275', 'ATP23', 'DISP3', 'HARs', 'B3GLCT', 'PLPP3', 'FGF7P3', 'EIPR1', 'MAIP1', 'VPS50', 'FLJ25363', 'GRAMD2B', '>', 'LRMDA', 'TEX45', 'CAVIN2', 'GAREM1', 'SLF1', 'TMEM131L', 'CFAP73', 'PLPP4', 'BAGE3', 'NOCT', 'PRKN', 'LNPK', 'DMAC1', 'NECTIN1', 'NSD3', 'LINC01600', 'SEM1', 'CCDC196']
#conversion to a pandas dataframe
deleted_genes = pd.DataFrame(deleted_genes, columns=['gene'])
#removing genes that are not found in the matrix from the gene set
geneset = geneset[~geneset.gene.isin(deleted_genes.gene)]
#       ---     NOW RUN CELL ABOVE AGAIN (LINES 30-32)     ---
#%% REMOVING GENES THAT HAVE ZERO EXPRESSION IN ALL THE CELLS
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
#If datasets are big and take a long time to load, here you can save the data. This will speed up loading times when working with the same dataset. 
#save the corrected datasets to csv genes
data_log10CPM.to_csv(r"path/example_name.csv")

#load saved (above) dataset
data_log10CPM = pd.read_csv(r"path/example_name.csv", index_col=0)

#%%ADDING AVERAGE EXPRESSION VALUES PER CELL
data_log10CPM['avg'] = data_log10CPM.mean(axis=1)

#%% IN CASE NEEDED TO NORMALIZE TO HOUSEKEEPING GENES, LOAD IN THE HOUSEKEEPING GENES DATASET
# load the HK values to calculate ratio between HK and gene set of interest  
HK_log10CPM = pd.read_csv(r"path/HK_geneset.csv", index_col=0)
average_HK = HK_log10CPM.mean(axis=1)
#%correct for HK ratio 
data_log10CPM['avg_correct'] = data_log10CPM['avg'] / average_HK

#%% ADDING METADATA TO THE DATASET
#add sample name, region and class information to the dataframe 
data_log10CPM['exp_component_name'] = list(samples_metadata_cleared['exp_component_name'])
data_log10CPM['region_label'] = list(samples_metadata_cleared['region_label'])
data_log10CPM['class_Label'] = list(samples_metadata_cleared['class_label'])

#%% SPLIT THE DATA BASED ON CORTICAL REGION 
MTG_avg = data_log10CPM.loc[data_log10CPM.region_label.values == 'MTG']
V1C_avg =  data_log10CPM.loc[data_log10CPM.region_label.values == 'V1C']
CgG_avg =  data_log10CPM.loc[data_log10CPM.region_label.values == 'CgG']
M1lm_avg = data_log10CPM.loc[data_log10CPM.region_label.values =='M1lm']
M1um_avg = data_log10CPM.loc[data_log10CPM.region_label.values =='M1ul']
S1ul_avg = data_log10CPM.loc[data_log10CPM.region_label.values =='S1ul']
S1lm_avg = data_log10CPM.loc[data_log10CPM.region_label.values =='S1lm']
A1C_avg =  data_log10CPM.loc[data_log10CPM.region_label.values =='A1C']
#from motor and somatosensory cortices merge the 2 subregions
M1_avg = pd.concat([M1lm_avg , M1um_avg], axis=0)
S1_avg = pd.concat([S1ul_avg, S1lm_avg], axis =0)

#%% SPLIT THE DATA INTO CELL CLASSES WITHIN SUB-REGIONS
#Glutamatergic
MTG_avg_Glu = MTG_avg.loc[(MTG_avg.class_Label.values == 'Glutamatergic')]
V1C_avg_Glu = V1C_avg.loc[(V1C_avg.class_Label.values == 'Glutamatergic')]
CgG_avg_Glu = CgG_avg.loc[(CgG_avg.class_Label.values == 'Glutamatergic')]
M1_avg_Glu = M1_avg.loc[(M1_avg.class_Label.values == 'Glutamatergic')]
S1_avg_Glu = S1_avg.loc[(S1_avg.class_Label.values == 'Glutamatergic')]
A1C_avg_Glu = A1C_avg.loc[(A1C_avg.class_Label.values == 'Glutamatergic')]

#GABAergic
MTG_avg_GABA = MTG_avg.loc[(MTG_avg.class_Label.values == 'GABAergic')]
V1C_avg_GABA = V1C_avg.loc[(V1C_avg.class_Label.values == 'GABAergic')]
CgG_avg_GABA = CgG_avg.loc[(CgG_avg.class_Label.values == 'GABAergic')]
M1_avg_GABA = M1_avg.loc[(M1_avg.class_Label.values == 'GABAergic')]
S1_avg_GABA = S1_avg.loc[(S1_avg.class_Label.values == 'GABAergic')]
A1C_avg_GABA = A1C_avg.loc[(A1C_avg.class_Label.values == 'GABAergic')]

#Non-neuronal
MTG_avg_non = MTG_avg.loc[(MTG_avg.class_Label.values == 'Non-neuronal')]
V1C_avg_non = V1C_avg.loc[(V1C_avg.class_Label.values == 'Non-neuronal')]
CgG_avg_non = CgG_avg.loc[(CgG_avg.class_Label.values == 'Non-neuronal')]
M1_avg_non = M1_avg.loc[(M1_avg.class_Label.values == 'Non-neuronal')]
S1_avg_non = S1_avg.loc[(S1_avg.class_Label.values == 'Non-neuronal')]
A1C_avg_non = A1C_avg.loc[(A1C_avg.class_Label.values == 'Non-neuronal')]

#%% PLOTTING SCRIPTS
"""
Sections below plot the data into violin plots with medians. 
All the axis labels, ranges and colors can be changed and adjusted by user. 

First section plots all the cells together per cortical region.
Following sections plot the data per cell class (Glutamatergic, GABAergic and non-Neuronal)

"""
#%% plotting expression values into violin plots, per sub-region for all cell classes. 
ax = sns.violinplot(data=[MTG_avg.avg, CgG_avg.avg, M1_avg.avg, V1C_avg.avg, S1_avg.avg, A1C_avg.avg], inner = 'quartile',scale="width", linewidth=4,
                   alpha=0.8)

ax.set_xticklabels(['MTG', "CgG","M1", "V1", "S1", "A1"])
ax.set_xlabel('brain region')
ax.set_ylabel("Mean Log10-CPM+1")
#ax.set_ylim(0, 1.5)
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

plt.gcf().set_size_inches(4,10)

plt.savefig("path/figure.eps", format='eps')
plt.savefig("path/figure.png", format ='png')

#%% plotting expression values into violin plots, per sub-region for glutamatergic cell-types.  
ax2 = sns.violinplot(data=[MTG_avg_Glu.avg, CgG_avg_Glu.avg, M1_avg_Glu.avg, V1C_avg_Glu.avg, S1_avg_Glu.avg, A1C_avg_Glu.avg], inner = 'quartile',scale="width", linewidth=4,
                   alpha=0.8, figsize=(10,16))
ax2.set_xticklabels(['MTG', "CgG","M1", "V1", "S1", "A1"])
ax2.set_xlabel('brain region')
ax2.set_ylabel("mean gene expression Log10-(CPM+1)")
#ax2.set_ylim(0,1.2)
for l in ax2.lines:
    l.set_linestyle('--')
    l.set_linewidth(0)
    l.set_color('red')
    l.set_alpha(0.8)
for l in ax2.lines[1::3]:
    l.set_linestyle('-')
    l.set_linewidth(3)
    l.set_color('black')
    l.set_alpha(0.8)
    
plt.gcf().set_size_inches(4,10)
    
plt.savefig("path/figure.eps", format='eps')
plt.savefig("path/figure.png", format ='png')

#%%plotting expression values into violin plots, per sub-region for GABAergic cell-types. 
ax3 = sns.violinplot(data=[MTG_avg_GABA.avg, CgG_avg_GABA.avg, M1_avg_GABA.avg, V1C_avg_GABA.avg, S1_avg_GABA.avg, A1C_avg_GABA.avg], inner = 'quartile',scale="width", linewidth=4,
                   alpha=0.8, figsize=(10,16))
ax3.set_xticklabels(['MTG', "CgG","M1", "V1", "S1", "A1"])
ax3.set_xlabel('brain region')
ax3.set_ylabel("Ratio EA/HK Mean Log10-(CPM+1)")
#ax3.set_ylim(0,4.7)
for l in ax3.lines:
    l.set_linestyle('--')
    l.set_linewidth(0)
    l.set_color('red')
    l.set_alpha(0.8)
for l in ax3.lines[1::3]:
    l.set_linestyle('-')
    l.set_linewidth(3)
    l.set_color('black')
    l.set_alpha(0.8)
 
plt.gcf().set_size_inches(4,10)
     
plt.savefig("path/figure.eps", format='eps')
plt.savefig("path/figure.png", format ='png')
    
#%%plotting expression values into violin plots, per sub-region for non-neuronal cell-types. 
ax4 = sns.violinplot(data=[MTG_avg_non.avg_correct, CgG_avg_non.avg_correct, M1_avg_non.avg_correct, V1C_avg_non.avg_correct, S1_avg_non.avg_correct, A1C_avg_non.avg_correct], inner = 'quartile',scale="width", linewidth=4,
                   alpha=0.8, figsize=(10,16))

ax4.set_xticklabels(['MTG', "CgG","M1", "V1", "S1", "A1"])
ax4.set_xlabel('brain region')
ax4.set_ylabel("Ratio EA/HK Mean Log10-(CPM+1)")
#ax4.set_ylim(0,4.7)
for l in ax4.lines:
    l.set_linestyle('--')
    l.set_linewidth(0)
    l.set_color('red')
    l.set_alpha(0.8)
for l in ax4.lines[1::3]:
    l.set_linestyle('-')
    l.set_linewidth(3)
    l.set_color('black')
    l.set_alpha(0.8)
    
#plt.gcf().set_size_inches(4,10)

plt.savefig("path/figure.eps", format='eps')
plt.savefig("npath/figure.png", format ='png')