

function [output] = Heatmap_plot_genes(Summ,cpm,var,data_genes,sign,lim, pos, neg)  
 % sign can be 'pos' for postively correlated genes and 'neg' for plotting negatively correlated genes   
% display options:
disp_genes_pos=pos;
disp_genes_neg=neg;
div=25;
% Make a table of gene expression of gene set and AP/TDL 
    cpm = cpm(ismember(cpm.sample_id,Summ.sample_id),:) ;
    Summ=Summ(ismember(Summ.sample_id,cpm.sample_id),:);
    % remove genes with zero expression in all cells
    indexes = [];
    
    for i = 3:numel(cpm.Properties.VariableNames)
        column = cpm{:,i};
        if all(column == 0)
            indexes = [indexes, i]; 
        end
    end
    cpm(:,indexes) = [];

    % select genes
    data_genes(isnan(data_genes.tStat),:)=[];
    data_genes=sortrows(data_genes,'tStat','descend');
    p_values=data_genes.p;
    R_values=data_genes.R2;
    genes=data_genes.Gene;
    genes = genes(ismember(genes,cpm.Properties.VariableNames)) ;
    select = cpm(:,['sample_id'; genes]) ;
    %%
       
    S=Summ(:,{'sample_id','t_type','NewLabel', var});        
    SummG = join(S,select,'Keys','sample_id') ;
    SummG=sortrows(SummG,var);
    SummG(isnan(SummG.(var)),:)=[];
%     
    SummG2=SummG(:,5:end);
    % get rid of genes with only 1 cell expression 
    % it is done in previous steps
%     indexes2=[];
%      for i = 1:size(SummG2,2)
%          column = table2array(SummG2(:,i));
%          if sum(column>0)==1
%              indexes2 = [indexes2, i]; 
%          end
%     end
%     SummG2(:, indexes2)=[];
    data=table2array(SummG2);
    
   
    %convert to log10(cpm+1)
    data=log10(data+1);
    data=data';
    yvalues = SummG2.Properties.VariableNames;
    xvalues=table2cell(SummG(:,'sample_id'));
    mid=floor(size(data,1)/2);
    
    for i=1:size(data,1)
       data2(i,:)=(data(i,:)-mean(data(i,:),'omitnan'))/std(data(i,:),'omitnan');
    end
  if strcmpi(sign,'pos')
   figure('Position',[312 129 1455 760]);
    % h = heatmap(xvalues(1:end),yvalues(1:disp_genes_pos),data(1:disp_genes_pos,1:end),...
   %       'ColorLimits',[0 lim]);
     h = heatmap(xvalues(1:end),yvalues(1:disp_genes_pos),data2(1:disp_genes_pos,1:end),...
          'ColorLimits',[-1 lim]);

    h.Colormap=jet;
    output=table(yvalues(1:disp_genes_pos)', p_values(1:disp_genes_pos),R_values(1:disp_genes_pos),...
        data2(1:disp_genes_pos,1:end));
    output.Properties.VariableNames = {'Gene', 'p-value', 'R2', 'z-scored gene expression'};
  elseif strcmpi(sign,'neg')
   figure('Position',[312 129 1455 380]);
    n_genes=size(data2,1)-(disp_genes_neg-1);
    h2 = heatmap(xvalues(1:end),yvalues(n_genes:end),data2(n_genes:end,1:end),...
        'ColorLimits',[-1 lim]);
    h2.Colormap=jet;
      output=table(yvalues(n_genes:end)', p_values(n_genes:end),R_values(n_genes:end),data2(n_genes:end,1:end));
      output.Properties.VariableNames = {'Gene', 'p-value', 'R2', 'z-scored gene expression'};

  end
 
    
