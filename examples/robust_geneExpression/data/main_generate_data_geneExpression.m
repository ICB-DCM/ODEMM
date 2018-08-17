clear all
close all
clc

distributions = {'norm','skew_norm','students_t','neg_binomial'};

n_cells = 100;
for iDist = 4%1:length(distributions)
    clear xi
    xi = [1,1,log10(20),0,log10(5),-1,0.3];
    if strcmp(distributions{iDist},'students_t')
        xi(8) = 3;
    end
    if strcmp(distributions{iDist},'skew_norm')
        xi(8) = 70;
    end
    D = generate_data_geneExpression('MA',distributions{iDist},xi,n_cells,0);
    plotODEMM(D)
    save(['data_geneExpression_MA_' distributions{iDist}],'D')
    
    D = generate_data_geneExpression('MA',distributions{iDist},xi,n_cells,0.1);
    save(['data_geneExpression_MA_' distributions{iDist} '_outlier'],'D')   
end

%% SSA
%D = generate_data_geneExpression('SSA',[],xi,n_cells);
