function D = generate_data_conversionProcess(mode,distribution,xi)

k1_1 = 10.^xi(1);
k1_2 = 10.^xi(2);
k2 = 10.^xi(3);
k3 = 10.^xi(4);
w = xi(5);

n_cells = 500;
t = [0,0.25,0.5,1];

D(1).t = t;
D(1).u = 1;
D(1).y = nan(1,numel(t),n_cells);

switch mode
    case 'MA'
        D(1).distribution = distribution;
        
        for s = 1:2
            if s == 1
                theta = [k1_1,k2,k3];
                ind_s = 1:round(n_cells*w);
                n_cells_s = length(ind_s);
            else
                theta = [k1_2,k2,k3];
                ind_s = (round(n_cells*w)+1):n_cells;
                n_cells_s = length(ind_s);
            end
            sol = simulate_conversionProcess_MA(t,theta);
            
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
                            mvtrnd(Sigma,nu,n_cells_s)+mu;
                        
                    case 'neg_binomial'
                        rho = 1 - m/C;
                        tau = (1 - rho)*m/rho;
                        
                        D(1).y(1,it,ind_s) = ...
                            nbinrnd(tau,rho,n_cells_s);
                end
            end
        end
        
    case 'SSA'
        warning off
        modelDefFileName = 'modelDef_conversionProcess';
        modelName = 'modelDef_conversionProcess';
        
        %method = 'MEC_2_LD_2_a';
        %System = genmexp(modelName,'modelDef_conversionProcess_foramici',method);
        %amiwrap('geneExp1D',[method,'_',modelName,'_syms'])
        
        eval(['System = ' modelDefFileName '([])']);
        System = completeSystem(System);
        System = completeSystemSSA(System);
        
        %options.mode = 'constant';
        
        %%
        theta_SSA_forinit = [0;k2;k3];
        % SSA runs for the initial conditions (without stimulation)
        [Xinits,Y_ssa,mX_ssa,mY_ssa,CX_ssa,CY_ssa] = simulate_SSA(System,[1,10,300],...
            theta_SSA_forinit,[],n_cells);
        
        theta_SSA = [k1_1;k2;k3];
        for k = 1:round(n_cells*w)
            System = modelDef_conversionProcess(Xinits(end,:,k)');
            System = completeSystem(System);
            System = completeSystemSSA(System);
            [X_ssa,Y_ssa,mX_ssa,mY_ssa,CX_ssa,CY_ssa] = simulate_SSA(System,t,theta_SSA,[],1);
            D(1).y(1,:,k) = permute(Y_ssa,[1,3,2]);
        end
        
        % Generation second population
        theta_SSA = [k1_2;k2;k3];
        for k = (round(n_cells*w)+1):n_cells
            System = modelDef_conversionProcess(Xinits(end,:,k)');
            System = completeSystem(System);
            System = completeSystemSSA(System);
            [X_ssa,Y_ssa,mX_ssa,mY_ssa,CX_ssa,CY_ssa] = simulate_SSA(System,t,theta_SSA,[],1);
            D(1).y(1,:,k) = permute(Y_ssa,[1,3,2]);
        end
        D(1).distribution = 'SSA';
end

D(1).n_dim = 1;
D(1).measurand = {'B'};
D(1).name = 'conversion reaction';
