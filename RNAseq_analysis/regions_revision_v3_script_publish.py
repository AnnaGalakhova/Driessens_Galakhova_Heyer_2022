# -*- coding: utf-8 -*-
"""
Created on Wed Mar 1 10:16:26 2023

This script is used to analyze and plot subsets of the AIBS RNA-seq database (https://portal.brain-map.org/atlases-and-data/rnaseq).
MULTIPLE CORTICAL AREAS - SMART-SEQ (2019) full matrix (exons+introns) is used to sub-select expression values
for specific genes of interest per region or cell classes. 
Script structure: 
24 - 38   selecting the data based on gene set of interest.
39 - 43   create datasets based on gene-set of interest (can be save and used later)    
44 - 61   Load in selected data and add values to metadata 
62 - 69   Group data by donor, allowing for plots of donordata seperate 
70 - 97  plot the data
98 - end statistical evalutation 

note. script on normaliztion of this dataset can be found in FINAL_RNA_seq_data_nromalization.py 

@author: Stan Driessens & Anna Galakhova 
"""
#load packages 
import pandas as pd 
import seaborn as sns 
import matplotlib.pyplot as plt
#read and import data from regions 
#import metadata
metadata = pd.read_csv(r'path/metadata.csv')
metadata = metadata.set_index('sample_name')
#clear the metadata
metadata_cleared = metadata[metadata['class_label'].notnull()]
#%% create gene-set data from normalized dataframe
geneset = pd.read_csv(r'path/geneset_IQ.txt', names=['gene'], header=None   , quotechar="'")
#in case there are duplicates present in the gene set, remove and update the gene set 
geneset = geneset.gene.drop_duplicates()
#convert geneset to a dataframe to improve compatibillity with pandas 
geneset = pd.DataFrame(geneset)
#%% load and filter data based on gene-set (ONLY NEED TO DO THIS ONCE, THEN SAVE PRE MADE DATA AS .CSV)
# selecting expression data of the genes present in the gene set of interest (geneset)
skiprows_set1 = metadata_cleared.loc[~metadata_cleared.sample_name.isin(metadata_cleared.sample_name)].index+1
#load in the data, NOTE. here we used pre-normalized log10 CPM data. 
data_log10CPM = pd.read_csv(r'path/matrix_log10_CPM.csv',
                  skiprows=skiprows_set1, usecols=geneset.gene)
#save the pre made data 
data_log10CPM.to_csv(r'path/your_data.csv')
#%%load in the pre-made normalized and pre-assigned gene-set data 
data = pd.read_csv(r'path/your_data.csv')
#create avg expression and add to metadata, align index on sample_name
data = data.set_index('sample_name')
#drop unnamed  0
data = data.drop('Unnamed: 0', axis=1)
data['avg'] = data.mean(axis=1)
#select only cpm values
data_cpm = data.avg 
#concatenate metadata to avg values from this geneset 
b = pd.concat([metadata_cleared, data_cpm], axis = 1 , sort = True)
#create new column for region values where M1 and S1 sub-regions are merged to new M1 and S1
b.loc[b['region_label'].str.contains('M1'), 'region']  = 'M1'
b.loc[b['region_label'].str.contains('MTG'), 'region']  = 'MTG'
b.loc[b['region_label'].str.contains('A1'), 'region']  = 'A1'
b.loc[b['region_label'].str.contains('CgG'), 'region']  = 'CgG'
b.loc[b['region_label'].str.contains('S1'), 'region']  = 'S1'
b.loc[b['region_label'].str.contains('V1'), 'region']  = 'V1'
#%% Sub-select data per donor (to plot median experssion per donor)
#here, select which cell class to plot (i.e. Glutamatergic, GABAergic or Non-neuronal)
data_plot  = b.loc[b['class_label'] == 'Glutamaterigc']
data_donor = b.loc[b['class_label'] == 'Glutamatergic']
#assign median expression per donor 
median_vals = data_donor.groupby(['region', 'external_donor_name_order', 'class_label'])['avg'].median()
median_df = pd.DataFrame(median_vals)
median_df.reset_index(inplace=True)
#%%plotting the data
#plot a violin for gene expression for all cells 
ax = sns.violinplot(data=data_plot, x='class_label', y='avg', legend=None, inner = None, 
                    order=['Glutamatergic', 'GABAergic', 'Non-neuronal'], color='lightgrey', linewidth = 0 )


ax = sns.pointplot(data = median_df, x='class_label', y='avg', order=['Glutamatergic', 'GABAergic', 'Non-neuronal'],
                   hue = 'external_donor_name_order', markers=['o','s','^'], color='k', join=False, linestyles='--', 
                   dodge=True,scale=1.3, errwidth=0)
    
ax.set_ylim((median_df.avg.mean()-0.25), median_df.avg.mean()+0.25)
ax.set_ylabel('Mean gene expression log10 (CPM+1)')

for l in ax.lines:
    l.set_linestyle('--')
    l.set_linewidth(0)
    l.set_color('red')
    l.set_alpha(0.5)
for l in ax.lines[1::3]:
    l.set_linestyle('--')
    l.set_linewidth(1)
    l.set_color('black')
    l.set_alpha(0.5)
plt.legend([],[], frameon=False)
#save the plot 
ax.set_ylim((median_df.avg.mean()-0.25), median_df.avg.mean()+0.25)
ax.set_ylabel('Mean gene expression log10 (CPM+1)')
plt.savefig(r'path/your_plot.eps')
#%%statistical evalutation 
#Friedman test on the separate donors
#split the group in separate donors 
#do Friedman test for geneset - class - regions 
# get the data per patient 
data_donor_1 = b.loc[b['external_donor_name_order'] == 1]
data_donor_2 = b.loc[b['external_donor_name_order'] == 2]
data_donor_3 = b.loc[b['external_donor_name_order'] == 3]
#get_class 
#hard code which class to use for statsitics: Glutamatergic, GABAergic or Non-neuronal 
data_donor_1_df = data_donor_1.loc[data_donor_1['class_label'] == 'Glutamatergic']
data_donor_2_df = data_donor_2.loc[data_donor_2['class_label'] == 'Glutamatergic']
data_donor_3_df = data_donor_3.loc[data_donor_3['class_label'] == 'Glutamatergic']

#%%do statistics, HERE SELECT WHICH DATAFRAME YOU WANT TO USE DONOR DATA OR COMBINED DATA
df = data_donor_1_df
#import staistical packages 
import scikit_posthocs as sp
import pingouin as pg
#perform Kruskal-Wallis and Dunn tests 
krus = pg.kruskal(data=df, dv='avg', between='region', detailed=True)
dunn = pg.pairwise_tests(dv='avg', between='region', data=df, padjust='holm', parametric=False)
#save the statistics as csv
krus.to_csv(r'path/kruskal_result.csv')
dunn.to_csv(r'path/dunn_result.csv')
#perform Friedman test 
fried = pg.friedman(data=df, dv='avg', between='region', within='donor')
dunn = pg.pairwise_tests(dv='avg', between='region', data=df, padjust='holm', parametric=False)
#save the statistics as csv
fried.to_csv(r'path/friedman_result.csv')
dunn.to_csv(r'path/dunn_fried_result.csv')

