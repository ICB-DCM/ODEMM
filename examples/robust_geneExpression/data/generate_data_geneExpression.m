function D = generate_data_geneExpression(mode,distribution,xi,n_cells,perc_outlier)

t = [0,5,10,30,60];
D(1).t = t;
D(1).u = 1;
D(1).y = NaN(1,numel(t),n_cells);

k1 = 10.^xi(1);
k2_1 = 10.^xi(2);
k2_2 = 10.^xi(3);
k3 = 10.^xi(4);
k4 = 10.^xi(5);
k5 = 10.^xi(6);
w = xi(7);

switch mode
    case 'MA'
        D(1).distribution = distribution;
        
        for s = 1:2
            if s == 1
                theta = [k1;k2_1;k3;k4;k5];
                ind_s = 1:round(n_cells*w);
                n_cells_s = length(ind_s);
            else
                theta = [k1;k2_2;k3;k4;k5];
                ind_s = (round(n_cells*w)+1):n_cells;
                n_cells_s = length(ind_s);
            end
            sol = simulate_geneExpression_MA(t,theta);
            
            for it = 1:numel(t)
                m = sol.y(it,1);
                C =  sol.y(it,2);
                switch distribution
                    case 'norm'
                        mu = m;
                        Sigma = C;
                        
                        D(1).y(1,it,ind_s) = ...
                            normrnd(mu,sqrt(Sigma),1,n_cells_s);
                        
                    case 'skew_norm'
                        delta = xi(end);
                        mu = m - sqrt(2/pi)*delta;
                        Sigma = C-(1-2/pi)*(delta*delta');
                        
                        D(1).y(1,it,ind_s) = ...
                            skewnormrnd(mu,Sigma,delta,n_cells_s);
                        
                    case 'students_t'
                        nu = xi(end);
                        mu = m;
                        Sigma = C;
                        
                        D(1).y(1,it,ind_s) = ...
                            studentstrnd(mu,Sigma,nu,n_cells_s);
                        
                    case 'neg_binomial'
                        rho = m/C;
                        tau = rho*m/(1-rho);
                        
                        D(1).y(1,it,ind_s) = ...
                            nbinrnd(tau,rho,n_cells_s,1);
                end
                if perc_outlier > 0
                    tmp = randperm(n_cells);
                    D(1).y(1,it,tmp(1:round(n_cells*perc_outlier))) = ...
                        randi(1800,round(n_cells*perc_outlier),1);
                end
            end
        end
        
        
    case 'SSA'
        modelDefFileName = 'modelDef_geneExpression';
        modelName = 'modelDef_geneExpression';
        
        method = 'MEC_2_LD_2_a';
        %System = genmexp(modelName,'modelDef_geneExp1D_foramici',method);
        %amiwrap('geneExp1D',[method,'_',modelName,'_syms'])
        
        eval(['System = ' modelDefFileName '([])']);
        System = completeSystem(System);
        System = completeSystemSSA(System);
        
        theta_SSA_forinit = [k1;0;k3;k4;k5];
        % SSA runs for the initial conditions (without stimulation)
        [Xinits,~,~,~,~,~] = simulate_SSA(System,[1,10,300],...
            theta_SSA_forinit,[],n_cells);
        
        theta_SSA = [k1;k2_1;k3;k4;k5];
        for k = 1:round(n_cells*w)
            System = modelDef_geneExperssion(Xinits(end,:,k)');
            System = completeSystem(System);
            System = completeSystemSSA(System);
            [~,Y_ssa,~,~,~,~] = simulate_SSA(System,t,theta_SSA,[],1);
            D(1).y(1,:,k) = permute(Y_ssa,[1,3,2]);
        end
        
        % Generation second population
        theta_SSA = [k1;k2_2;k3;k4;k5];
        for k = (round(n_cells*w)+1):n_cells
            System = modelDef_geneExpression(Xinits(end,:,k)');
            System = completeSystem(System);
            System = completeSystemSSA(System);
            [~,Y_ssa,~,~,~,~] = simulate_SSA(System,t,theta_SSA,[],1);
            D(1).y(1,:,k) = permute(Y_ssa,[1,3,2]);
        end
        D(1).distribution = 'SSA';
end

D(1).n_dim = 1;
D(1).measurand = {'A'};
D(1).name = 'Gene expression';
D(1).xi_true = xi;

% figure
% theta_SSA = [k1;k2;k3;k4;k5];
% sol = simulate_geneExp1D(D.t,theta_SSA);
% plot(sol.t,sol.y(:,1),'r-'); hold on;
% plot(D.t, squeeze(mean(D(1).y(1,:,1:30),3)),'b--');
% theta_SSA = [k1;k2_2;k3;k4;k5];
% sol = simulate_geneExp1D(D.t,theta_SSA);
% plot(sol.t,sol.y(:,1),'r-'); hold on;
% plot(D.t, squeeze(mean(D(1).y(1,:,31:end),3)),'b--');