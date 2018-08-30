clear all
close all
clc

% Generate all ODEMMs for all distribution assumptions and varying number
% of timepoints
%generateODEMMs_robust_geneExpression()

%%
cellNumbers = [50,100,500,1000,5000];
timePoints{1} = [0,1,2];
timePoints{2} = [0,1,2,5];
timePoints{3} = [0,1,2,4,5];
timePoints{4} = [0,1,2,3,4,5];

distributions = {'norm','skew_norm','students_t','neg_binomial'};

for iDist = 1:4
    for iCells = 1:length(cellNumbers)
        for iTime = 1:length(timePoints)
            for iDataset = 1:5
                for jDist = 1:4
                    if ~(jDist == 4 && iDist <4) % neg binomial only for neg binomial data
                        clear parameters M D options conditions
                        load([distributions{iDist} '/data_geneExpression_simstudy_MA_' ...
                            distributions{iDist} ...
                            '_' num2str(cellNumbers(iCells)) 'cells_' ...
                            num2str(length(timePoints{iTime})) 'tps_' num2str(iDataset)],'D')
                        
                        if ~isequal(length(timePoints{iTime}),4)
                            M.name = [distributions{jDist} '_geneExpression_2_' num2str(length(timePoints{iTime})) 'tps' ];
                        else
                            M.name = [distributions{jDist} '_geneExpression_2' ];
                        end
                        eval(['ODEMM_' M.name]);
                        if strcmp(distributions{jDist},'skew_norm')
                            parameters.min(end) = -2e2;
                            parameters.max(end) = 2e2;
                        end
                        if strcmp(distributions{jDist},'students_t')
                            parameters.max(end) = 8;
                        end
                        [conditions,D] = collectConditions(D,M);

                        options.MS = PestoOptions();
                        options.MS.n_starts = 30;
                        options.MS.localOptimizerOptions.Display = 'silent';
                        
                        parameters.guess = getParameterGuesses(parameters,@(xi) ...
                            logLikelihood_extend(xi,M,D,options,conditions),...
                            options.MS.n_starts, parameters.min,parameters.max);
                        warning off
                        parameters = getMultiStarts(parameters,@(xi) ...
                            logLikelihood_extend(xi,M,D,options,conditions),options.MS);
                        save(['./results/' distributions{iDist} ...
                            '/results_geneExpression_simstudy_MA_' ...
                            distributions{iDist} ...
                            '_' num2str(cellNumbers(iCells)) 'cells_' ...
                            num2str(length(timePoints{iTime})) 'tps_' num2str(iDataset) ...
                            '' distributions{jDist}])
                    end
                end
            end
        end
    end
end
