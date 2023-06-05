# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 10:43:50 2023
This script is used to analysze and plot subsets of the AIBS RNA-seq database (https://portal.brain-map.org/atlases-and-data/rnaseq).
MTG - SMART-SEQ (2019) seperate matrices (exons and introns) are used to sub-select expression values
for specific genes of interest per transcriptomic type within MTG. 

Script structure: 
17 - 33  selecting metadata
35 - 60  load and obtain data for genesets of interest, including CPM+1 log10 transformation and saving the data
61 - 84 sturcturing data for further plotting
85 - 108  plotting
109 - end Statistics
@authors: Stan Driessens & Anna Galakhova 
"""

#load packages 
import pandas as pd 
import seaborn as sns 
import numpy as np 
import matplotlib.pyplot as plt


#read and import data for MTG RNA-seq data  
#import metadata
metadata = pd.read_csv(r'path/human_MTG_2018-06-14_samples-columns.csv')
metadata_cleared = metadata.set_index('sample_name')
#%%make data set per geneset
gene_selection = pd.read_csv(r'path/geneset.txt',
                             names=['gene'], header=None   , quotechar="'")

usecols = list(samples_metadata_cleared.sample_name)
usecols.append('Unnamed: 0')

#load in reads from exon and intron 
exons = pd.read_csv(r"path/human_MTG_2018-06-14_exon-matrix.csv", 
                        usecols=usecols, index_col = 0)
introns = pd.read_csv(r"path/human_MTG_2018-06-14_exon-matrix.csv", 
                        usecols=usecols, index_col = 0)
#add exon + intron reads 
data = exons + introns
#transpose the data to make genes columns 
data_T = data.T
#create a column with total reads per cell
data_T['total_reads'] = data_T.sum(axis=1)
#divide all reads by total reads 
data_T_c = data_T.div(data_T['total_reads'], axis=0)
#multiply by 1e6 to get 'counts per million'
data_T_cpm = data_T_c*1e6
#add 1 to prevent inf values when log transforming
data_T_cpm1 = data_T_cpm+1
#transorm 10 log 
data_log10cpm = np.log10(data_T_cpm1)
#check if the correct genes are in the dataframe 
selection_i = gene_metadata[gene_metadata.gene.isin(gene_selection.gene)]
selection = list(selection_i.entrez_id)
#subselect the data to only contain genes of interest
output = data_log10cpm[selection]
#%% save the pre made data 
output.to_csv(r'path/my_data.csv')
#%%load pre made data set
data = pd.read_csv(r'path/my_data.csv', index_col = 'Unnamed: 0')
#get an average expresssion value for each cell 
data_cpm = data.mean(axis=1)
data_cpm = data_cpm.rename('avg')
#add average data to metadata creating one big long format Sort = True to allign to correct index 
final_data = pd.concat([metadata_cleared, data_cpm], axis = 1, sort=True)
#only select layer 2/3 from the data
final_data_23 = final_data.loc[(final_data['brain_subregion'] == 'L2') | (final_data['brain_subregion'] == 'L3')]
#select the t-types of interest
final_data_23_types = final_data_23[(final_data_23['cluster'].str.contains('FREM3')) | (final_data_23['cluster'].str.contains('GLP2'))
                  | (final_data_23['cluster'].str.contains('LTK')) | (final_data_23['cluster'].str.contains('CARM1P1'))
                  | (final_data_23['cluster'].str.contains('COL22'))]
#subselect frem into deep and superficial 
final_data_23_types['cluster_new']=list(map(lambda x: x.split()[3], final_data_23_types.cluster))
final_data_23_types.loc[(final_data_23_types.cluster_new.values == 'FREM3') & (final_data_23_types.brain_subregion.values == 'L2'), 'cluster_new'] = 'L2 FREM3'
final_data_23_types.loc[(final_data_23_types.cluster_new.values == 'FREM3') & (final_data_23_types.brain_subregion.values == 'L3'), 'cluster_new'] = 'L3 FREM3'


#%%create median values for stripplot purpose (per donor)
median_vals = final_data_23_types.groupby(['donor', 'cluster_new'])['avg'].median()
median_df = pd.DataFrame(median_vals)
median_df.reset_index(inplace=True)

#%%PLOT
sns.set_palette(sns.color_palette('pastel'))

ax = sns.violinplot(data=final_data_23_types, x='cluster_new', y='avg', legend=None, inner = 'quartile',
                    order=['LTK', 'L2 FREM3', 'L3 FREM3', 'GLP2R', 'CARM1P1', 'COL22A1'], color='lightgrey')
ax = sns.pointplot(data = median_df, x='cluster_new', y='avg', order=['LTK', 'L2 FREM3', 'L3 FREM3' , 'GLP2R', 'CARM1P1', 'COL22A1'],
                   hue = 'donor', markers=['o','s','^', 'p'], color='k', join=False, linestyles='-', 
                   dodge=True,scale=1.3, errwidth=0)
ax.set_ylim(0, 1.75)
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
ax.set_ylim(0.6, 1.3)
ax.set_ylabel('Mean gene expression log10 (CPM+1)')
#%%save the plot
plt.savefig(r'path/my_plot.eps')
#%%STATISTICS
#%%do statisics for every donor seperately 
donors = np.unique(final_data_23_types.donor)
H16_24_010 = final_data_23_types[final_data_23['donor'] == 'H16.24.010']
H200_1023 = final_data_23_types[final_data_23_types['donor'] == 'H200.1023']
H200_1025 = final_data_23_types[final_data_23_types['donor'] == 'H200.1025']
H200_1030 = final_data_23_types[final_data_23_types['donor'] == 'H200.1030']
#%%do statistics, HERE SELECT WHICH DATAFRAME YOU WANT TO USE DONOR DATA OR COMBINED DATA
#import statistical packages 
import scikit_posthocs as sp
import pingouin as pg
#do Kruskal-Wallis and Dunn for single cell data
krus = pg.kruskal(data=df, dv='avg', between='cluster_new', detailed=True)
dunn = pg.pairwise_tests(dv='avg', between='cluster_new', data=df, padjust='holm', parametric=False)
#save the statistics as csv
krus.to_csv(r'path/kruskal_result.csv')
dunn.to_csv(r'path/dunn_result.csv')
#perforam Friedmann test for data per donor
fried = pg.friedman(data=df, dv='avg', between='region', within='donor')
dunn = pg.pairwise_tests(dv='avg', between='region', data=df, padjust='holm', parametric=False)
#save the statistics as csv
fried.to_csv(r'path/friedman_result.csv')
dunn.to_csv(r'path/dunn_fried_result.csv')
