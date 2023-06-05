basedir='/Users/natalia/Documents/Manuscripts/2023 Driessens Genes of intelligence/Final files after revision/Scripts /Scripts Analysis Patch-seq data';
cpm=readtable(fullfile(basedir,'gene_cmp_all.csv'));
load data.mat
%%
[p, tbl, comp, NMIQR, X] = dotplot_genes(SummH2,cpm,{'SCN2B'}, 'AvgExpLog10',[]);
%%
[p,N,R2,tStat]=showcorrGeneset(SummH2,cpm,{'CACNA2D3'},'AvgExpLog10','upstrokeFB0',[],[],'Upstroke')
%%
[p,N,R2,tStat]=showcorrGeneset(SummH2,cpm,{'CADPS2'},'AvgExpLog10','upstroke21to40',[],[],'Upstroke')
%%
[p,N,R2,tStat]=showcorrGeneset(SummHmorph2,cpm,{'HTR2A'},'AvgExpLog10','TDL',[],[],'TDL, um')
%%
genes=cpm.Properties.VariableNames(7:end)';
%%
 AP=MostSigGenes2(SummH2,'All_cells',cpm,genes,'upstroke21to40'); 
 %writetable(AP,fullfile(basedir,'MSG_AP_first_all_genes.csv'));

  %%

 AP_IQ=MostSigGenes(SummH2,'All_cells',cpm,gsets(1),'upstroke21to40'); 
  %writetable(stats,fullfile(savedir,'MostSigGenes/MSG_IQ_AP_first_all_cells.csv'));

 AP_EA=MostSigGenes(SummH2,'All_cells',cpm,gsets(4),'upstroke21to40') ;
   %writetable(stats,fullfile(savedir,'MostSigGenes/MSG_EA_AP_first_all_cells.csv'));
 
 AP_HAR=MostSigGenes(SummH2,'All_cells',cpm,gsets(5),'upstroke21to40') ;
  %writetable(stats,fullfile(savedir,'MostSigGenes/MSG_HAR_AP_first_all_cells.csv'));
%%
test1=AP_IQ1(ismember(AP_IQ1,AP_EA1),:);
test2=AP_EA1(ismember(AP_EA1,AP_H1),:);
test3=AP_IQ1(ismember(AP_IQ1,AP_H1),:);