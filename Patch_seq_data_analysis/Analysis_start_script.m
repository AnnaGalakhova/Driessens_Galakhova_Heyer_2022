% last update 2023-04-18 N Goriounova
clear all
close all

%% Analysis script

%% load data 

% This files contains morphology. physiology, gene expression, cell type data 
% cpm values file (counts per million) for genes that are expressed in 5%
% of all cells
cpm_all_data=readtable('cpm_all_data.csv');
% the last 3 columns of cpm_all_data are: {'TDL_um'}
% {'AP_rise_speed_mV…'}    {'t_type'}. Delete this to have a table of only
% cpm values
cpm=cpm_all_data(:, 1:size(cpm_all_data,2)-3);
data=cpm_all_data(:, {'sample_id','TDL_um', 'AP_rise_speed_mVms', 't_type', 'NewLabel', 'donor_label'});
% genesets.mat contains the 3 genesets: IQ, EA and HAR genes
load genesets.mat
%
%% 
% for i=1:size(cpm_all_data,1)
%     idx=strcmpi(cpm_all_data.sample_id(i), SummH2.sample_id);
%     cpm_all_data.donor_label(i)=SummH2.donor_label(idx);
% end

%% Figure 2C plots
ylimit=[0.6 1.4];

[p,tbl,comp,NMIQR, X_IQ] = dotplot_Geneset(data,cpm,IQ_genes, 'AvgExpLog10',ylimit); 

ylimit=[0.8 1.8];
[p,tbl,comp,NMIQR, X_EA] = dotplot_Geneset(data,cpm,EA_genes, 'AvgExpLog10',ylimit); 

ylimit=[0.6 1.7];
[p,tbl,comp,NMIQR, X_HAR] = dotplot_Geneset(data,cpm,HAR_genes, 'AvgExpLog10',ylimit); 
%print(fullfile(savedir,'Figure 2/Genes_H_in_types.eps'), '-depsc','-r300','-vector');


%% Figure 3A (TDL)
[stats_TDL_types, NMIQR, p, X_TDL] = dotplottype_NG(data, 'TDL_um', [],'Total dendritic length, mm');
%print(fullfile(savedir,'Figure 3/TDL_types.eps'), '-depsc','-r300','-painters');

 %% Figure Supplementary 2: dotplots per donor
% gene expression per donor
ylimit=[0.6 1.4];
[p,tbl,comp,NMIQR, X_IQ] = dotplot_Geneset_per_donor(data,cpm,IQ_genes, 'AvgExpLog10',ylimit) ;
%print(fullfile(savedir,'Figure Suppl Donor/Genes_S_in_types_donor.eps'), '-depsc','-r300','-vector');
ylimit=[0.8 1.8];
[p,tbl,comp,NMIQR, X_EA] = dotplot_Geneset_per_donor(data,cpm,EA_genes, 'AvgExpLog10',ylimit) 

%print(fullfile(savedir,'Figure Suppl Donor/Genes_L_in_types_donor.eps'), '-depsc','-r300','-vector');
ylimit=[0.8 1.6];

[p,tbl,comp,NMIQR, X_H] = dotplot_Geneset_per_donor(data,cpm,HAR_genes, 'AvgExpLog10',ylimit) 
%print(fullfile(savedir,'Figure Suppl Donor/Genes_H_in_types_donor.eps'), '-depsc','-r300','-vector');

%
ylimit=[0 20000]
[stats_TDL_types, NMIQR, p, X_TDL] = dotplottype_per_donor(data, 'TDL_um', ylimit,'Total dendritic length, mm');
%print(fullfile(savedir,'Figure Suppl Donor/TDL_types_donor.eps'), '-depsc','-r300','-painters');
%
ylimit=[0 600]
[stats_AP_types,NMIQR, p, X_AP] = dotplottype_per_donor(data, 'AP_rise_speed_mVms', ylimit,'First AP rise speed , mV/ms');
%print(fullfile(savedir,'Figure Suppl Donor/AP_first_types_donor.eps'), '-depsc','-r300','-painters');

%% Figure 4C
ylimit=[0 600]
[stats_AP_types,NMIQR, p, X_AP] = dotplottype_NG(data, 'AP_rise_speed_mVms', ylimit,'First AP rise speed , mV/ms');
%print(fullfile(savedir,'Figure 4/AP_first_types.eps'), '-depsc','-r300','-painters');
 


%% identify most significant single genes to use for GSEA 


 TDL_IQ=MostSigGenes(data,cpm,IQ_genes,'TDL_um') ;
 
 TDL_EA=MostSigGenes(data,cpm,EA_genes,'TDL_um') ;
 
 TDL_HAR=MostSigGenes(data,cpm,HAR_genes,'TDL_um') ;

%%
 AP_IQ=MostSigGenes(data,cpm,IQ_genes,'AP_rise_speed_mVms'); 

 AP_EA=MostSigGenes(data,cpm,EA_genes,'AP_rise_speed_mVms') ;
 
 AP_HAR=MostSigGenes(data,cpm,HAR_genes,'AP_rise_speed_mVms') ;
 
%% Load all genes 

% remove genes that with NAN values
TDL_IQ(isnan(TDL_IQ.p),:)=[]; 
TDL_EA(isnan(TDL_EA.p),:)=[];
TDL_HAR(isnan(TDL_HAR.p),:)=[];


AP_IQ(isnan(AP_IQ.p),:)=[]; 
AP_EA(isnan(AP_EA.p),:)=[];
AP_HAR(isnan(AP_HAR.p),:)=[];
%% Get the significant values after multiple comparisons correction

% run correction for multiple testing using Benjamini & Hochberg procedure
% and FDR of 0.2?
FDR=0.05;
[h, crit_p_IQ_AP, adj_ci_cvrg, adj_p]=fdr_bh(AP_IQ.p,FDR,'pdep','yes');

[h, crit_p_EA_AP, adj_ci_cvrg, adj_p]=fdr_bh(AP_EA.p,FDR,'pdep','yes');

[h, crit_p_HAR_AP, adj_ci_cvrg, adj_p]=fdr_bh(AP_HAR.p,FDR,'pdep','yes');

[h, crit_p_IQ_TDL, adj_ci_cvrg, adj_p]=fdr_bh(TDL_IQ.p,FDR,'pdep','yes');

[h, crit_p_EA_TDL, adj_ci_cvrg, adj_p]=fdr_bh(TDL_EA.p,FDR,'pdep','yes');

[h, crit_p_HAR_TDL, adj_ci_cvrg, adj_p]=fdr_bh(TDL_HAR.p,FDR,'pdep','yes');

%% select the positively and negatively correlating genes
AP_IQ_genes=AP_IQ(AP_IQ.p<=crit_p_IQ_AP,:);
AP_EA_genes=AP_EA(AP_EA.p<=crit_p_EA_AP,:);
AP_HAR_genes=AP_HAR(AP_HAR.p<=crit_p_HAR_AP,:);
TDL_IQ_genes=TDL_IQ(TDL_IQ.p<=crit_p_IQ_TDL,:);
TDL_EA_genes=TDL_EA(TDL_EA.p<=crit_p_EA_TDL,:);
TDL_HAR_genes=TDL_HAR(TDL_HAR.p<=crit_p_HAR_TDL,:);

%% make the heatmaps for Figure 3C,D,E of positively and negatively correlated genes with TDL
% 
[output_TDL_IQ_pos]=Heatmap_plot_genes(data,cpm,'TDL_um', TDL_IQ,'pos',3, 33, 4);
%print(fullfile(savedir,'Heatmaps/IQ_TDL_pos.eps'), '-depsc','-r300','-vector');
[output_TDL_IQ_neg]=Heatmap_plot_genes(data,cpm,'TDL_um', TDL_IQ,'neg',3, 33, 4);
%print(fullfile(savedir,'Heatmaps/IQ_TDL_neg.eps'), '-depsc','-r300','-vector');

[output_TDL_EA_pos]=Heatmap_plot_genes(data,cpm,'TDL_um', TDL_EA,'pos',3, 56, 12);
%print(fullfile(savedir,'Heatmaps/EA_TDL_pos.eps'), '-depsc','-r300','-vector');

[output_TDL_EA_neg]=Heatmap_plot_genes(data,cpm,'TDL_um', TDL_EA,'neg',3, 56, 12);
%print(fullfile(savedir,'Heatmaps/EA_TDL_neg.eps'), '-depsc','-r300','-vector');

[output_TDL_HAR_pos]=Heatmap_plot_genes(data,cpm,'TDL_um', TDL_HAR,'pos',3, 72, 13);
%print(fullfile(savedir,'Heatmaps/HAR_TDL_pos.eps'), '-depsc','-r300','-vector');

[output_TDL_HAR_neg]=Heatmap_plot_genes(data,cpm,'TDL_um', TDL_HAR,'neg',3, 72, 13);
%print(fullfile(savedir,'Heatmaps/HAR_TDL_neg.eps'), '-depsc','-r300','-vector');
%% %% make the heatmaps for Figure 4E, F, G of positively and negatively correlated genes with AP
% AP correlations
[output_AP_IQ_pos]=Heatmap_plot_genes(data,cpm,'AP_rise_speed_mVms', AP_IQ,'pos',3, 23, 3);
%print(fullfile(savedir,'Heatmaps/IQ_AP_pos.eps'), '-depsc','-r300','-vector');

[output_AP_IQ_neg]=Heatmap_plot_genes(data,cpm,'AP_rise_speed_mVms', AP_IQ,'neg',3, 23, 3);
%print(fullfile(savedir,'Heatmaps/IQ_AP_neg.eps'), '-depsc','-r300','-vector');

[output_AP_EA_pos]=Heatmap_plot_genes(data,cpm,'AP_rise_speed_mVms', AP_EA,'pos',3, 39, 10);
%print(fullfile(savedir,'Heatmaps/EA_AP_pos.eps'), '-depsc','-r300','-vector');

[output_AP_EA_neg]=Heatmap_plot_genes(data,cpm,'AP_rise_speed_mVms', AP_EA,'neg',3, 39, 10);
%print(fullfile(savedir,'Heatmaps/EA_AP_neg.eps'), '-depsc','-r300','-vector');

[output_AP_HAR_pos]=Heatmap_plot_genes(data,cpm,'AP_rise_speed_mVms', AP_HAR,'pos',3, 71, 10);
%print(fullfile(savedir,'Heatmaps/HAR_AP_pos.eps'), '-depsc','-r300','-vector');

[output_AP_HAR_neg]=Heatmap_plot_genes(data,cpm,'AP_rise_speed_mVms', AP_HAR,'neg',3, 71, 10);
%print(fullfile(savedir,'Heatmaps/HAR_AP_neg.eps'), '-depsc','-r300','-vector');


%% Display average expression in correlated genes per type
%% Figure 3C
gset=TDL_IQ_genes;

genes=table2cell(gset(gset.tStat>0,1));
[p, tbl, comp, NMIQR, X] = dotplot_selected_genes(data,cpm, genes, []);

genes=table2cell(gset(gset.tStat<0,1));
[p, tbl, comp, NMIQR, X] = dotplot_selected_genes(data,cpm, genes, []);
%% Figure 3D
gset=TDL_EA_genes;
genes=table2cell(gset(gset.tStat>0,1));
[p, tbl, comp, NMIQR, X] = dotplot_selected_genes(data,cpm, genes, []);
%print(fullfile(savedir,'Figure 3/Pos_types_L.eps'), '-depsc','-r300','-painters');
genes=table2cell(gset(gset.tStat<0,1));
[p, tbl, comp, NMIQR, X] = dotplot_selected_genes(data,cpm, genes, []);
%print(fullfile(savedir,'Figure 3/Neg_types_L.eps'), '-depsc','-r300','-painters');
%% Figure 3E
gset=TDL_HAR_genes;
genes=table2cell(gset(gset.tStat>0,1));
[p, tbl, comp, NMIQR, X] = dotplot_selected_genes(data,cpm, genes, []);
%print(fullfile(savedir,'Figure 3/Pos_types_H.eps'), '-depsc','-r300','-painters');

genes=table2cell(gset(gset.tStat<0,1));
[p, tbl, comp, NMIQR, X] = dotplot_selected_genes(data,cpm, genes, []);
%print(fullfile(savedir,'Figure 3/Neg_types_H.eps'), '-depsc','-r300','-painters');

%% Figure 4E
gset=AP_IQ_genes;
genes=table2cell(gset(gset.tStat>0,1));
[p, tbl, comp, NMIQR, X] = dotplot_selected_genes(data,cpm, genes, []);
%print(fullfile(savedir,'Figure 4/Pos_types_S.eps'), '-depsc','-r300','-painters');

genes=table2cell(gset(gset.tStat<0,1));
[p, tbl, comp, NMIQR, X] = dotplot_selected_genes(data,cpm, genes, []);
%print(fullfile(savedir,'Figure 4/Neg_types_S.eps'), '-depsc','-r300','-painters');
%% Figure 4F
gset=AP_EA_genes;
genes=table2cell(gset(gset.tStat>0,1));
[p, tbl, comp, NMIQR, X] = dotplot_selected_genes(data,cpm, genes, []);
%print(fullfile(savedir,'Figure 4/Pos_types_L.eps'), '-depsc','-r300','-painters');
genes=table2cell(gset(gset.tStat<0,1));
[p, tbl, comp, NMIQR, X] = dotplot_selected_genes(data,cpm, genes, []);
%print(fullfile(savedir,'Figure 4/Neg_types_L.eps'), '-depsc','-r300','-painters');
%% Figure 4G
gset=AP_HAR_genes;
genes=table2cell(gset(gset.tStat>0,1));
[p, tbl, comp, NMIQR, X] = dotplot_selected_genes(data,cpm,genes, []);
%print(fullfile(savedir,'Figure 4/Pos_types_H.eps'), '-depsc','-r300','-painters');

genes=table2cell(gset(gset.tStat<0,1));
[p, tbl, comp, NMIQR, X] = dotplot_selected_genes(data,cpm, genes, []);
%print(fullfile(savedir,'Figure 4/Neg_types_H.eps'), '-depsc','-r300','-painters');


%% Select genes that correlated significantly with AP and TDL (after correction)

% find overlap sig genes
AP_overlap_SL=AP_IQ_genes(ismember(AP_IQ_genes, AP_EA_genes),:);
AP_overlap_SH=AP_IQ_genes(ismember(AP_IQ_genes, AP_HAR_genes),:);
AP_overlap_LH=AP_EA_genes(ismember(AP_EA_genes, AP_HAR_genes),:);
AP_overlap_all=AP_IQ_genes(ismember(AP_IQ_genes,AP_EA_genes)&ismember(AP_IQ_genes,AP_HAR_genes),:);
[C,ia] = unique(AP_overlap_all.Gene);
AP_overlap = AP_overlap_all(ia,:);


TDL_overlap_SL=TDL_IQ_genes(ismember(TDL_IQ_genes, TDL_EA_genes),:);
TDL_overlap_SH=TDL_IQ_genes(ismember(TDL_IQ_genes, TDL_HAR_genes),:);
TDL_overlap_LH=TDL_EA_genes(ismember(TDL_EA_genes, TDL_HAR_genes),:);
TDL_overlap_all=TDL_IQ_genes(ismember(TDL_IQ_genes,TDL_EA_genes)&ismember(TDL_IQ_genes,TDL_HAR_genes),:);
[C,ia2] = unique(TDL_overlap_all.Gene);
TDL_overlap = TDL_overlap_all(ia2,:);
%% for figure 4; find overlap between AP and TDL correlated genes together

overlap=outerjoin(TDL_overlap, AP_overlap, 'Keys', 'Gene');

overlap_S=outerjoin(TDL_IQ_genes,AP_IQ_genes, 'Keys', 'Gene'); 
overlap_L=outerjoin(TDL_EA_genes, AP_EA_genes, 'Keys', 'Gene'); 
overlap_H=outerjoin(TDL_HAR_genes, AP_HAR_genes, 'Keys', 'Gene'); 
overlap_SL=outerjoin(TDL_overlap_SL, AP_overlap_SL, 'Keys', 'Gene'); 
overlap_SH=outerjoin(TDL_overlap_SH, AP_overlap_SH, 'Keys', 'Gene');
overlap_LH=outerjoin(TDL_overlap_LH, AP_overlap_LH, 'Keys', 'Gene');
overlap_all=outerjoin(TDL_overlap_all, AP_overlap_all, 'Keys', 'Gene');

% find overlap in the complete gene sets
AP_overlap_SL2=AP_IQ(ismember(AP_IQ, AP_EA),:);
AP_overlap_SH2=AP_IQ(ismember(AP_IQ, AP_HAR),:);
AP_overlap_LH2=AP_EA(ismember(AP_EA, AP_HAR),:);
AP_overlap_all2=AP_IQ(ismember(AP_IQ,AP_EA)&ismember(AP_IQ,AP_HAR),:);
[C,ia] = unique(AP_overlap_all2.Gene);
AP_overlap2 = AP_overlap_all2(ia,:);


TDL_overlap_SL2=TDL_IQ(ismember(TDL_IQ, TDL_EA),:);
TDL_overlap_SH2=TDL_IQ(ismember(TDL_IQ, TDL_HAR),:);
TDL_overlap_LH2=TDL_EA(ismember(TDL_EA, TDL_HAR),:);
TDL_overlap_all2=TDL_IQ(ismember(TDL_IQ,TDL_EA)&ismember(TDL_IQ,TDL_HAR),:);
[C,ia2] = unique(TDL_overlap_all2.Gene);
TDL_overlap = TDL_overlap_all2(ia2,:);
%% %%

%% Make pie charts Figure 5

clearvars d1 d2 d3 d4 X

d2=size(overlap_SH,1);
d1=size(overlap_S,1)-d2;

X=[d1 d2];
figure;
donut(X,'ringwidth',0.5);
colormap(lines(2))
title('IQ genes (Savage)significantly correlated with AP and TDL')
%print(fullfile(savedir,'Figure 5/Donut_Savage_sig.eps'), '-depsc','-r300','-vector');

 %
clearvars d1 d2 d3 d4 X
d2=size(overlap_LH,1);
d1=size(overlap_L,1)-d2;


X=[d1 d2];
figure;

donut(X,'ringwidth',0.5);
colormap(lines(2));
title('EA genes (Lee) significantly correlated with AP and TDL')
%print(fullfile(savedir,'Figure 5/Donut_Lee_sig.eps'), '-depsc','-r300','-vector');



%% Venn chart of overlap Figure 5
clear vars xlim ylim

figure('Position', [100 200 500 450]); 
% for all genes:
venn([size(overlap_H,1) size(overlap_S,1) size(overlap_L,1)], ...
    [size(overlap_SH,1) size(overlap_LH,1) size(overlap_SL,1) size(overlap_all,1) ],...
    'FaceColor',{'y','r','b'})

xlim([-8 12])
ylim([-8 12])

%print(fullfile(savedir,'Figure 5/Venn_Sig_AP_TDL.eps'), '-depsc','-r300','-vector');
