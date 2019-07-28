clear all
close all
clc

method = 'MEC_2_LD_2_a';
modelnames = {'conversionReaction','diffProteinExpression',...
    'twoStageGeneExpression','diffProteinExpression'};

for m = [1,3,4]
    modelDefFileName = ['modelDef_model' num2str(m) '_' modelnames{m} '_foramici'];
    modelName = ['model' num2str(m) '_' modelnames{m}];
    System_MM2 = genSimFile(modelName,modelDefFileName,method);    
end