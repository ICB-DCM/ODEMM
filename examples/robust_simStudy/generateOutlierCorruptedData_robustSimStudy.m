% This script generates the outlier-corrupted data for the simulation study.

clear all
close all
clc

% Models included in the manuscript number 1, 3 and 4
modelnames = {'conversionReaction',... % 1D data
    'diffProteinExpression',... % 2D data
    'twoStageGeneExpression',... % 1D data
    'diffProteinExpression'}; % 1D data

load_simStudy_settings

% outlierPerc = 0.02;
% outlierstr = 'outlier2_zeros';
 
% outlierPerc = 0.05;
% outlierstr = 'outlier5_dublets';

outlierPerc = 0.1;
outlierstr = 'outlier10_unif';

% If folder does not exist, generate it
if ~(exist('dataOutlier')==7)
    mkdir dataOutlier
end

for m = [1,3,4]
    for it = 1:3
        t = tps{it};
        for ic = 1:4
            for set = 1:3
                load(['./data/data_model' num2str(m) '_' modelnames{m} '_' ...
                    num2str(n_cells(ic)) ...
                    'cells_' num2str(length(t)) 'tps_' num2str(set) 'paramsetD'],'D')
                nOutlier = round(n_cells(ic)*outlierPerc);
                indOutlier = randperm(n_cells(ic));
                indOutlier = indOutlier(1:nOutlier);
                switch outlierstr
                    case 'outlier2_zeros'
                        for iit = 1:length(t)
                            D.y(:,iit,indOutlier,:) = 0;
                        end
                    case 'outlier5_dublets'
                        for iit = 1:length(t)
                            D.y(:,iit,indOutlier,:) = 2*D.y(:,iit,indOutlier,:);
                        end
                    case 'outlier10_unif'
                        for iit = 1:length(t)
                            interv = max(squeeze(D.y(:,iit,:,:)))-min(squeeze(D.y(:,iit,:,:)));
                            lb = max(0,min(squeeze(D.y(:,iit,:,:)))-interv*0.25);
                            ub = max(squeeze(D.y(:,iit,:,:)))+interv*0.25;
                            D.y(:,iit,indOutlier,:) = round(rand(length(indOutlier),1)*(lb-ub)+ub);
                        end
                end
                save(['./dataOutlier/data_model' num2str(m) '_' ...
                    modelnames{m} '_' num2str(n_cells(ic)) ...
                    'cells_' num2str(length(t)) 'tps_' num2str(set) ...
                    'paramset_' outlierstr],'D')
            end
        end
    end
end
