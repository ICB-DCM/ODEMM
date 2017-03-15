% Generation of all models accounting for subpopulation differences

close all
clear all

load('./project/data/data_Lys_1D2D_sepscaled');

parameters.max = [6*ones(9,1);... %kinetic
    4*ones(3,1);...%correlation cv
    4*ones(8,1);...%scaling
    1*ones(6,1); % noise
    1]; % weight
parameters.min = [-6*ones(9,1);... %kinetic
    -4*ones(3,1);...%scaling
    -4*ones(8,1);...%scaling
    -6*ones(6,1); % noise
    0]; % weight
    

differences = {'log_{10}({k_1}_{s2})','log_{10}({k_2}_{s2})','log_{10}({k_4}_{s2})','log_{10}({k_5}_{s2})',...
     'log_{10}({k_3TrkA_0}_{s2})','log_{10}({s_{P1}Erk}_{s2})'};


comb = nan(2^6,6);
i = 1;
for i1 = 0:1
   for i2 = 0:1
    for i3 = 0:1
        for i4 = 0:1
            for i5 = 0:1
                for i6 = 0:1
                    comb(i,:) = [i1,i2,i3,i4,i5,i6];
                    i = i+1;
                end
            end
        end
    end
   end 
end
% [0 0 0 1 1 0]
% [0 0 0 1 0 0]

for icomb = [8     9    25    32    34    39    41    48    50    64]
parameters.name = {'log_{10}(k_1)','log_{10}(k_2)','log_{10}(k_4)','log_{10}(k_5)',...
    'log_{10}(k_3TrkA_0)','log_{10}(s_{P1}Erk)',...
    'log_{10}(cvTrkA)','log_{10}(cvErk)','log_{10}(corr)'...
    'log_{10}(b_{P1})','log_{10}(s_{P2})','log_{10}(b_{P2})',...
    'log_{10}(s_{P3})','log_{10}(b_{P3})','log_{10}(s_{P4})','log_{10}(b_{P4})',...
    'log_{10}(sT_3/k3)','log_{10}(bT_3)',...
    'log_{10}(s_{E4}/s_{P1})','log_{10}(b_{E4})',...
    'log_{10}(\epsilon_P_{1})','log_{10}(\epsilon_P_{2})','log_{10}(\epsilon_P_{3})',...
    'log_{10}(\epsilon_P_{4})','log_{10}(\epsilon_T_{3})','log_{10}(\epsilon_E_{4})',...
    'w'};
M.n_subpop = 2; %number of subpopulations
M.model = @(T,theta,u)simulate_SigmaPoints_sPsET_loglog_corr(T,theta,u(1)); 
% 
% P E T PP PE PT EE ET TT
for s = 1:M.n_subpop
for e=1:2
    M.mean_ind{s,e} = [1]; %index of mean in output
    M.var_ind{s,e} = [4]; %index of variances in output (empty if RREs used)
end
e = 3;
M.mean_ind{s,e} = [3,1]; %TrkA/pErk
M.var_ind{s,e} = [9,6,4];
e = 4;
M.mean_ind{s,e} = [2,1]; %Erk/pErk
M.var_ind{s,e} = [7,5,4]; 

end

M.sim_type = 'HO'; 

options.replicates = false;
options.write_parameters = true;
options.measurement_noise = true;
options.noise_model = 'multiplicative';
options.penalize_scaling = false;

inddiff = 28;
str_symtheta = 'M.sym.theta = [';
    for id = 1:6
             
     if comb(icomb,id) 
         str_symtheta = [str_symtheta 'log(10)*(u(g+1)*xi(' num2str(id) ')'];     
         str_symtheta = [str_symtheta '+(1-u(g+1))*xi(' num2str(inddiff) ')'];
         parameters.name{end+1} = differences{id};
         parameters.max(end+1) = 6;
         parameters.min(end+1) = -6;
         inddiff = inddiff+1;
     else
         str_symtheta = [str_symtheta 'log(10)*(xi(' num2str(id) ')'];  
     end
     str_symtheta = [str_symtheta ');'];
    end
str_symtheta = [str_symtheta 'log(log(10.^(2*(xi(7)))+1));log(10)*(xi(9));log(log(10.^(2*(xi(8)))+1))];'];
parameters.number = length(parameters.name);

% definition of symbolic parameters required for generateODEMM
xi = sym(zeros(parameters.number,1));
for i = 1:parameters.number
    xi(i) = sym(['xi_' num2str(i,'%d')]);
end
u = sym(zeros(2,1)); % 1 for definition of stimulation, 1 for subpop,
for i = 1:7
   u(i) = sym(['u_' num2str(i,'%d')]);
end
g = size(D(e).u,1); % number of input for simulation of inidividual subpopulation

for e=[1:4]       
    M.u{1,e} = [0]; % subpop
    M.u{2,e} = [1];
end

r=1; 
e=1;
M.sym.scaling{r,e} = sym('1');
M.sym.offset{r,e} = 10.^(xi(10));
M.sym.sigma_noise{e} = 10.^(xi(21));

e=2;
M.sym.scaling{r,e} = 10.^(xi(11));
M.sym.offset{r,e} = 10.^(xi(12));
M.sym.sigma_noise{e} = 10.^(xi(22));

e=3;
M.sym.scaling{r,e} = [10.^xi(17);10.^xi(13)];
M.sym.offset{r,e} = [10.^xi(18);10.^xi(14)];
M.sym.sigma_noise{e} = [10.^xi(25);10.^xi(23)];

e=4;
M.sym.scaling{r,e} = [10.^xi(19);10.^xi(15)];
M.sym.offset{r,e} = [10.^xi(20);10.^xi(16)];
M.sym.sigma_noise{e} = [10.^xi(26);10.^xi(24)];

for e = 1:4
    for s = 1:2
        M.distribution{s,e} = 'logn_mean'; % 'norm', 'logn_median'
    end
    % initialize weights
    M.sym.w{1,e} = ones(size(D(e).t(:)))*xi(27);
    M.sym.w{2,e} = ones(size(D(e).t(:)))*(1-xi(27));
end   
eval(str_symtheta);
M.name = ['NGF_subpop_comb' num2str(icomb) ]; 


[conditions,D] = collectConditions(D,M);
M = generateODEMM(D,M,parameters,conditions,options);

end


