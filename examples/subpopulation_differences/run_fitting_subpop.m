function parameters = run_fitting_subpop(icomb)


load('./data/data_PDL');
options.MS = PestoOptions();
options.MS.n_starts = 100;
options.MS.localOptimizer = 'fmincon';
options.MS.localOptimizerOptions = optimset('GradObj','on',...
    'display','iter','TolFun',1e-10, 'TolX',1e-14, 'MaxIter', 1000,...
    'algorithm','interior-point');
options.MS.comp_type = 'sequential';
options.MS.mode = 'text';
options.simulate_musigma = true;

try
    eval(['ODEMM_NGF_subpop_comb' num2str(icomb)]);
catch
    error(['ODEMM_NGF_subpop_comb' num2str(icomb) 'not found. Run generate_subpop_files.m to generate the model.'])
end

options.MS.foldername = ['results_' M.name];
[conditions,D] = collectConditions(D,M);

parameters.guess = getParameterGuesses(parameters,@(xi) ...
    logLikelihood(xi,M,D,options,conditions),...
    options.MS.n_starts,parameters.min,parameters.max);

parameters = getMultiStarts(parameters,@(xi) ...
    logLikelihood(xi,M,D,options,conditions),options.MS);

