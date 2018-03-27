% In this script the data for poly-D-lysine is loaded and saved such that 
% it can be used within the ODEMM framework. 

clear all;
close all;
clc;

ExpC1 = load_NGF_pErk_kin('PDL');
s1 = getScalingFactors('log',ExpC1);

ExpC2 = load_NGF_pErk_dr('PDL');
s2 = getScalingFactors('log',ExpC2);

ExpC3 = load_NGF_dose_response_TrkApErk('PDL');
s3 = getScalingFactors('log',ExpC3);

ExpC4 = load_NGF_dose_response_ErkpErk('PDL');
s4 = getScalingFactors('log',ExpC4);


%% Fill struct D needed for ODEMM likelihood
% NOTE: For the model simulation the NGF concentration is considered in a
% different unit (multiplication with 1/20)

ExpC = ExpC1;
D(1).name = 'NGF kinetic - PDL';
D(1).n_dim = 1;
D(1).measurand = 'pErk levels';
D(1).t = [0,1,5,15,30,60,120];
D(1).u = 1;
D(1).y = nan(length(D(1).u),length(D(1).t),6000);
D(1).Ey = nan(length(D(1).u),length(D(1).t),10);
for r = 1:3
    D(1).replicate{r} = nan(length(D(1).u),length(D(1).t),6000);
end
maxn0 = 0;
maxr = 0;
for k = 1:length(D(1).t)
    n0 = 0;
    nr = 0;
    for r = 1:length(ExpC(k).replicate)
        nx = size(ExpC(k).replicate(r).data{1},1);
        D(1).y(1,k,n0+1:n0+nx) = s1(r)*ExpC(k).replicate(r).data{1}(1:nx,1)';
        D(1).replicate{r}(1,k,1:nx) = s1(r)*ExpC(k).replicate(r).data{1}(1:nx,1)';
        D(1).Ey(1,k,r) = s1(r)*mean(ExpC(k).replicate(r).data{1}(1:nx,1)');
        n0 = n0 + nx;
        nr = nr +1;
    end
    
    if maxn0 < n0
        maxn0 = n0;
    end
    if maxr < nr
        maxr = nr;
    end
end
temp = D(1).y;
D(1).y = temp(:,:,1:maxn0);
D(1).Ey = D(1).Ey(:,:,1:maxr);
clear ExpC
%%
ExpC = ExpC2;
D(2).name = 'NGF dose response - PDL';
D(2).measurand = 'relative pErk conc.';
D(2).n_dim = 1;
D(2).t = [60];
D(2).u = [0,0.008,0.04,0.2,1,5,25];
D(2).y = nan(size(D(2).u,2),length(D(2).t),10000);
D(2).Ey = nan(size(D(2).u,2),length(D(2).t),10);
for r = 1:4
    D(2).replicate{r} = nan(length(D(2).u),length(D(2).t),6000);
end

maxn0=0;
maxr = 0;
for d = 1:size(D(2).u,2)
    n0 = 0;
    nr = 0;
    for r = 1:length(ExpC(d).replicate)
        nx = size(ExpC(d).replicate(r).data{1},1);
        D(2).y(d,1,n0+1:n0+nx) = s2(r)*ExpC(d).replicate(r).data{1}(1:nx,1)';
        D(2).replicate{r}(d,1,1:nx) = s2(r)*ExpC(d).replicate(r).data{1}(1:nx,1)';
        D(2).Ey(d,1,r) = s2(r)*mean(ExpC(d).replicate(r).data{1}(1:nx,1)');
        n0 = n0 + nx;
      nr = nr +1;
    end
    if maxn0 < n0
        maxn0 = n0;
    end
    if maxr < nr
        maxr = nr;
    end
end
temp = D(2).y;
D(2).y = temp(:,:,1:maxn0);
D(2).Ey = D(2).Ey(:,:,1:maxr);

%%
ExpC = ExpC3;
D(3).name = 'NGF dose response - TrkApErk - PDL';
D(3).n_dim = 2;
D(3).measurand = {'rel. TrkA conc.', 'rel. pErk conc.'};
D(3).t = 60;
D(3).u = [0,0.008,0.04,0.2,1,5];
D(3).y = nan(length(D(3).u),1,10000,2);
D(3).Ey = nan(length(D(3).u),1,10,2);
for r = 1:4
    D(3).replicate{r} = nan(length(D(3).u),length(D(3).t),6000,2);
end

maxn0 = 0;
maxr = 0;
for d = 1:length(D(3).u)
    n0 = 0;
    for r = 1:length(ExpC(d).replicate)
        nx = size(ExpC(d).replicate(r).data{1},1);
        D(3).y(d,1,n0+1:n0+nx,:) = s3(r)*ExpC(d).replicate(r).data{1}(1:nx,:);
        D(3).replicate{r}(d,1,1:nx,:) = s3(r)*ExpC(d).replicate(r).data{1}(1:nx,:);
        D(3).Ey(d,1,r,:) = s3(r)*mean(ExpC(d).replicate(r).data{1}(1:nx,:));
        n0 = n0 + nx;
    end
    if n0>maxn0
      maxn0 = n0;
    end
    if r > maxr
        maxr = r;
    end
end
D(3).y = D(3).y(:,:,1:maxn0,:);
D(3).Ey = D(3).Ey(:,1,1:maxr,:);

%%
ExpC = ExpC4;
D(4).name = 'NGF dose response - ErkpErk - PDL';
D(4).n_dim = 2;
D(4).measurand = {'rel. Erk conc.', 'rel. pErk conc.'};
D(4).t = 60;
D(4).u = [0,0.008,0.04,0.2,1,5];
D(4).y = nan(length(D(4).u),1,10000,2);
D(4).Ey = nan(length(D(4).u),1,10,2);
for r = 1:4
    D(4).replicate{r} = nan(length(D(4).u),length(D(4).t),6000,2);
end

maxn0 = 0;
maxr = 0;
for d = 1:length(D(4).u)
    n0 = 0;
    for r = 1:length(ExpC(d).replicate)
        nx = size(ExpC(d).replicate(r).data{1},1);
        D(4).y(d,1,n0+1:n0+nx,:) = s4(r)*ExpC(d).replicate(r).data{1}(1:nx,:);
        D(4).replicate{r}(d,1,1:nx,:) = s4(r)*ExpC(d).replicate(r).data{1}(1:nx,:);
        D(4).Ey(d,1,r,:) = s4(r)*mean(ExpC(d).replicate(r).data{1}(1:nx,:));
        n0 = n0 + nx;
    end

    if n0>maxn0
      maxn0 = n0;
    end
    if r > maxr
        maxr = r;
    end
end
D(4).y = D(4).y(:,:,1:maxn0,:);
D(4).Ey = D(4).Ey(:,1,1:maxr,:);

for e=1:4
    D(e).y(D(e).y<0) = NaN;
    for r = 1:numel(D(e).replicate)
        D(e).replicate{r}(D(e).replicate{r}<0) = NaN;
    end
end

%% save data
save data_PDL D