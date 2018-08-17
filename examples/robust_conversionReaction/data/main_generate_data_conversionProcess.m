clear all
close all
clc

distributions = {'norm','skew_norm','students_t','neg_binomial'};

xi = [-0.1,0.1,-0.45,-0.2,0.3];

for iDist = 1:length(distributions)
    if strcmp(distributions{iDist},'students_t')
        xi(7) = 300;
    end
    if strcmp(distributions{iDist},'skew_norm')
        xi(7) = 22;
    end
    D = generate_data_conversionProcess('MA',distributions{iDist},xi);
    plotODEMM(D)
    save(['data_conversionProcess_MA_' distributions{iDist}],'D')
end

%%
D = generate_data_conversionProcess('SSA',[],xi);
