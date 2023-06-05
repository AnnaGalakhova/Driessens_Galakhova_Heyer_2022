function [stats] = MostSigGenes(data,cpm,genes,var)  
%     if ismember(var, {'upstroke1to20', 'upstroke21to40','upstrokeFB0'})
%         Summ = Summ(ismember(Summ.collaborator_label,'ZZ_Missing'),:) ;
%     end
    genes = genes(ismember(genes,cpm.Properties.VariableNames)) ;
    select = cpm(ismember(cpm.sample_id,data.sample_id),['sample_id'; genes]) ;
    select{:,2:width(select)} = log10(select{:,2:width(select)}+1) ;
    cpm = cpm(ismember(cpm.sample_id,data.sample_id),:) ;
    data=data(ismember(data.sample_id,cpm.sample_id),:);

    SummG = join(data,select,'Keys','sample_id') ;
    SummG(isnan(SummG.(var)),:)=[];
    for i = 1:length(genes)   
        mdl = fitlm(SummG.(genes{i}), SummG.(var));
        stats(i).Gene = genes{i} ;
        stats(i).R2 = mdl.Rsquared.Ordinary;
        stats(i).p = mdl.Coefficients.pValue(2);
        stats(i).tStat = mdl.Coefficients.tStat(2);
       
    end
    stats = sortrows(struct2table(stats),3) ;
    
end



