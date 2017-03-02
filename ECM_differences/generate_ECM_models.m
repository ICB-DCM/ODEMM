function [] = generate_ECM_models()
% This function generates the 128 hierarchical models accounting for all
% possible differences between the extracellular scaffolds.
% The underlying model for each scaffolds assumes inter- and
% intra-subpopulation variability of TrkA activity and cell-to-cell
% variability of relative Erk1/2 levels.

load('./project/data/data_matrices_1D2D');

parameters.max = [6*ones(7,1);... %kinetic
    4*ones(4,1); ... % corr
    6*ones(11,1);...%scaling
    1*ones(6,1); % noise
    0]; % weight
parameters.min = [-6*ones(7,1);... %kinetic
    -4*ones(4,1); ... % corr
    -6*ones(11,1);...%scaling
    -3*ones(6,1); % noise
    -4]; %weight

differences = {'log_{10}(\kappa_{k_1})','log_{10}(\kappa_{k_2})','log_{10}(\kappa_{k_4})','log_{10}(\kappa_{k_5})',...
    'log_{10}(\kappa_{k_3[TrkA]_0})','log_{10}(\kappa_{(c_P[Erk]_{0})','log_{10}(\kappa_{w})'};

comb = nan(2^7,7);
i = 1;
for i1 = 0:1
    for i2 = 0:1
        for i3 = 0:1
            for i4 = 0:1
                for i5 = 0:1
                    for i6 = 0:1
                        for i7 = 0:1
                            comb(i,:) = [i1,i2,i3,i4,i5,i6,i7];
                            i = i+1;
                        end
                    end
                end
            end
        end
    end
end
%%
for icomb = 1:128
    
    parameters.name = {'log_{10}(k_1)','log_{10}(k_2)','log_{10}(k_4)','log_{10}(k_{5})',...% parameter for simulation
        'log_{10}(\beta_{k_3[TrkA]_{0,1}})','log_{10}(\beta_{k_3TrkA_{0,2}})',...
        'log_{10}(c_P^{1,2}[Erk]_{0,1})',...
        'log_{10}(\sigma_{T,1})','log_{10}(\sigma_{T,2})','log_{10}(\sigma_{E})','log_{10}(\sigma_{TE})'...
        'log_{10}(c_P^{1,2})','log_{10}(c_T^{5,6}/k_3)','log_{10}(o_T^{5,6})',... % scaling and offset
        'log_{10}(c_E^{7,8}/c_P^{5,6})','log_{10}(o_E^{7,8})',...
        'log_{10}(c_P^{3,4})','log_{10}(o_P^{3,4})',...
        'log_{10}(c_P^{5,6})','log_{10}(o_P^{5,6})',...
        'log_{10}(c_P^{7,8})','log_{10}(o_P^{7,8})',...
        'log_{10}(\sigma_{P,noise}^{1,2})','log_{10}(\sigma_{P,noise}^{3,4})',... %measurement noise
        'log_{10}(\sigma_{T,noise}^{5,6})','log_{10}(\sigma_{P,noise}^{5,6})',...
        'log_{10}(\sigma_{E,noise}^{7,8})','log_{10}(\sigma_{P,noise}^{7,8})',...
        'log_{10}(w)'};
    
    M.n_subpop = 2; %number of subpopulations
    M.model = @(T,theta,u)simulate_SigmaPoints_sPsET_loglog_corr(T,theta,u(1));
    %
    % P E T PP PE PT EE ET TT
    for s = 1:M.n_subpop
        for e=1:4
            M.mean_ind{s,e} = [1]; %index of mean in output
            M.var_ind{s,e} = [4]; %index of variances in output (empty if RREs used)
        end
        for e = 5:6
            M.mean_ind{s,e} = [3,1];
            M.var_ind{s,e} = [9,6,4];
        end
        for e=7:8
            M.mean_ind{s,e} = [2,1];
            M.var_ind{s,e} = [7,5,4];
        end
    end
    
    M.sim_type = 'HO';
    
    
    options.replicates = false;
    options.write_parameters = true;
    options.measurement_noise = true;
    options.noise_model = 'multiplicative';
    options.penalize_scaling = false;
    
    inddiff = 30;
    str_symtheta = 'M.sym.theta = [';
    for id = 1:4
        str_symtheta = [str_symtheta 'log(10)*(xi(' num2str(id) ')'];
        if comb(icomb,id)
            str_symtheta = [str_symtheta '+u(g+2)*xi(' num2str(inddiff) ')'];
            parameters.name{end+1} = differences{id};
            parameters.max(end+1) = 3;
            parameters.min(end+1) = -3;
            inddiff = inddiff+1;
        end
        str_symtheta = [str_symtheta ');'];
    end
    id = 5;
    str_symtheta = [str_symtheta 'log(10)*(u(g+1)*xi(5)+ (1-u(g+1))*xi(6)'];
    if comb(icomb,id)
        str_symtheta = [str_symtheta '+u(g+2)*xi(' num2str(inddiff) ')'];
        parameters.name{end+1} = differences{id};
        parameters.max(end+1) = 3;
        parameters.min(end+1) = -3;
        inddiff = inddiff+1;
    end
    str_symtheta = [str_symtheta ');'];
    id = 6;
    str_symtheta = [str_symtheta 'log(10)*(xi(7)'];
    if comb(icomb,id)
        str_symtheta = [str_symtheta '+u(g+2)*xi(' num2str(inddiff) ')'];
        parameters.name{end+1} = differences{id};
        parameters.max(end+1) = 3;
        parameters.min(end+1) = -3;
        inddiff = inddiff+1;
    end
    str_symtheta = [str_symtheta ');'];
    
    if comb(icomb,7)
        parameters.number = length(parameters.name) +1;
    else
        parameters.number = length(parameters.name);
    end
    xi = sym(zeros(parameters.number,1));
    for i = 1:parameters.number
        xi(i) = sym(['xi_' num2str(i,'%d')]);
    end
    u = sym(zeros(3,1)); % 1 for definition of stimulation, 1 for subpop, 1 for experiment
    for i = 1:7
        u(i) = sym(['u_' num2str(i,'%d')]);
    end
    g = size(D(e).u,1);
    
    if comb(icomb,7)
        for e = 1:8
            % initialize weights
            M.sym.w{1,e} = ones(size(D(e).t(:)))*(10.^(xi(29)+u(g+2)*xi(inddiff)));
            M.sym.w{2,e} = ones(size(D(e).t(:)))*(1-(10.^(xi(29)+u(g+2)*xi(inddiff))));
        end
        parameters.name{end+1} = differences{7};
        parameters.max(end+1) = 3;
        parameters.min(end+1) = -3;
        parameters.constraints.A = [zeros(28,1);1;zeros(inddiff-29-1,1);1];
        parameters.constraints.b = 0;
    else
        for e = 1:8
            M.sym.w{1,e} = ones(size(D(e).t(:)))*(10.^xi(29));
            M.sym.w{2,e} = ones(size(D(e).t(:)))*(1-10.^xi(29));
        end
    end
    str_symtheta = [str_symtheta 'log(log(10.^(2*(u(g+1)*xi(8) + (1-u(g+1))*xi(9)))+1));log(10)*xi(11);log(log(10.^(2*(xi(10)))+1))];'];
    
    parameters.number = length(parameters.name);
    
    for e = 1:8
        for s = 1:2
            M.distribution{s,e} = 'logn_mean'; % 'norm', 'logn_median'
        end
    end
    
    for e=[1,3,5,7]
        M.u{1,e} = [0;0]; % subpop, condition
        M.u{2,e} = [1;0];
    end
    for e=[2,4,6,8]
        M.u{1,e} = [0;1];
        M.u{2,e} = [1;1];
    end
    noise_ind = 23;
    r=1;
    for e=1:2
        M.sym.scaling{r,e} = sym('1');
        M.sym.offset{r,e} = 10.^(xi(12));
        M.sym.sigma_noise{e} = 10.^(xi(noise_ind));
    end
    noise_ind = noise_ind+1;
    for e=3:4
        M.sym.scaling{r,e} = 10.^(xi(17));
        M.sym.offset{r,e} = 10.^(xi(18));
        M.sym.sigma_noise{e} = 10.^(xi(noise_ind));
    end
    noise_ind = noise_ind+1;
    
    for e = 5:6
        M.sym.scaling{r,e} = [10.^(xi(13));10.^(xi(19))];
        M.sym.offset{r,e} = [10.^(xi(14));...
            10.^xi(20)];
        M.sym.sigma_noise{e} = [10.^(xi(noise_ind)); 10.^(xi(noise_ind+1))];
    end
    noise_ind = noise_ind+2;
    
    for e = 7:8
        M.sym.scaling{r,e} = [10.^xi(15);10.^(xi(21))];
        M.sym.offset{r,e} = [10.^xi(16);...
            10.^xi(22)];
        M.sym.sigma_noise{e} = [10.^(xi(noise_ind)); 10.^(xi(noise_ind+1))];
    end
    eval(str_symtheta);
    M.name = ['NGF_ECM_T_comb' num2str(icomb)];
    [conditions,D] = collectConditions(D,M);
    M = generateODEMM(D,M,parameters,conditions,options);
    eval(['ODEMM_NGF_ECM_T_comb' num2str(icomb)])
    options.simulate_musigma = true;
end
