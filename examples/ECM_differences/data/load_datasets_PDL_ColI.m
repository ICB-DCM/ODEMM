% Script for loading the datasets of extracellular scaffolds poly-D-lysin 
% (PDL) and Collagen I (ColI).

clear all;
close all;
clc;

ExpC1 = load_NGF_pErk_kin('PDL','ColI');
s1 = getScalingFactors('log',ExpC1);

ExpC2 = load_NGF_pErk_dr('PDL','ColI');
s2 = getScalingFactors('log',ExpC2);

ExpC3 = load_NGF_dose_response_TrkApErk('PDL','ColI');
s3 = getScalingFactors('log',ExpC3);
 
ExpC4 = load_NGF_dose_response_ErkpErk('PDL','ColI');
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
for r = 1:length(ExpC(1).replicate)
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
        D(1).Ey(1,k,r) = s1(r)*mean(ExpC(k).replicate(r).data{1}(1:nx,1)');
        D(1).replicate{r}(1,k,1:nx) = s1(r)*ExpC(k).replicate(r).data{1}(1:nx,1)';
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

D(2).name = 'NGF kinetic - ColI';
D(2).n_dim = 1;
D(2).measurand = 'pErk levels';
D(2).t = [0,1,5,15,30,60,120];
D(2).u = 1;
D(2).y = nan(length(D(2).u),length(D(2).t),6000);
D(2).Ey = nan(length(D(2).u),length(D(2).t),10);
for r = 1:length(ExpC(1).replicate)
	D(2).replicate{r} = nan(length(D(2).u),length(D(2).t),6000);
end

maxn0 = 0;
maxr = 0;
for k = 1:length(D(2).t)
    n0 = 0;
    nr = 0;
    for r = 1:length(ExpC(k).replicate)
        nx = size(ExpC(k).replicate(r).data{2},1);
        D(2).y(1,k,n0+1:n0+nx) = s1(r)*ExpC(k).replicate(r).data{2}(1:nx,1)';
        D(2).Ey(1,k,r) = s1(r)*mean(ExpC(k).replicate(r).data{2}(1:nx,1)');
      	D(2).replicate{r}(1,k,1:nx) = s1(r)*ExpC(k).replicate(r).data{2}(1:nx,1)';
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

clear ExpC

ExpC = ExpC2;
D(3).name = 'NGF dose response - PDL';
D(3).measurand = 'relative pErk conc.';
D(3).n_dim = 1;
D(3).t = [60];
D(3).u = [0,0.008,0.04,0.2,1,5,25];
D(3).y = nan(size(D(3).u,2),length(D(3).t),10000);
D(3).Ey = nan(size(D(3).u,2),length(D(3).t),10);
for r = 1:length(ExpC(1).replicate)
	D(3).replicate{r} = nan(size(D(3).u,2),length(D(3).t),6000);
end

maxn0=0;
for d = 1:size(D(3).u,2)
    n0 = 0;
    nr = 0;
    for r = 1:4
        nx = size(ExpC(d).replicate(r).data{1},1);
        D(3).y(d,1,n0+1:n0+nx) = s2(r)*ExpC(d).replicate(r).data{1}(1:nx,1)';
        D(3).Ey(d,1,r) = s2(r)*mean(ExpC(d).replicate(r).data{1}(1:nx,1)');
        D(3).replicate{r}(d,1,1:nx) = s2(r)*ExpC(d).replicate(r).data{1}(1:nx,1)';
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
temp = D(3).y;
D(3).y = temp(:,:,1:maxn0);
D(3).Ey = D(3).Ey(:,:,1:maxr);

D(4).name = 'NGF dose response - ColI';
D(4).measurand = 'relative pErk conc.';
D(4).n_dim = 1;
D(4).t = [60];
D(4).u = [0,0.008,0.04,0.2,1,5,25];
D(4).y = nan(size(D(4).u,2),length(D(4).t),10000);
D(4).Ey = nan(size(D(4).u,2),length(D(4).t),10);
for r = 1:length(ExpC(1).replicate)
	D(4).replicate{r} = nan(size(D(4).u,2),length(D(4).t),6000);
end

maxn0=0;
for d = 1:size(D(4).u,2)
    n0 = 0;
    nr = 0;
    for r = 1:4
        nx = size(ExpC(d).replicate(r).data{2},1);
        D(4).y(d,1,n0+1:n0+nx) = s2(r)*ExpC(d).replicate(r).data{2}(1:nx,1)';
        D(4).Ey(d,1,r) = s2(r)*mean(ExpC(d).replicate(r).data{2}(1:nx,1)');
        D(4).replicate{r}(d,1,1:nx) = s2(r)*ExpC(d).replicate(r).data{2}(1:nx,1)';
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
temp = D(4).y;
D(4).y = temp(:,:,1:maxn0);
D(4).Ey = D(4).Ey(:,:,1:maxr);

ExpC = ExpC3;
D(5).name = 'NGF dose response - TrkApErk - PDL';
D(5).n_dim = 2;
D(5).measurand = {'rel. TrkA conc.', 'rel. pErk conc.'};
D(5).t = 60;
D(5).u = [0,0.008,0.04,0.2,1,5];
D(5).y = nan(length(D(5).u),1,10000,2);
D(5).Ey = nan(length(D(5).u),1,10,2);
for r = 1:length(ExpC(1).replicate)
	D(5).replicate{r} = nan(length(D(5).u),1,10000,2);
end

n0max = 0;
rmax = 0;
for d = 1:length(D(5).u)
    n0 = 0;
    for r = 1:4
        nx = size(ExpC(d).replicate(r).data{1},1);
        D(5).y(d,1,n0+1:n0+nx,:) = s3(r)*ExpC(d).replicate(r).data{1}(1:nx,:);
        D(5).Ey(d,1,r,:) = s3(r)*mean(ExpC(d).replicate(r).data{1}(1:nx,:));
        D(5).replicate{r}(d,1,1:nx,:) = s3(r)*ExpC(d).replicate(r).data{1}(1:nx,:);
        n0 = n0 + nx;
    end
        if n0>n0max
      n0max = n0;
    end
    if r > rmax
        rmax = r;
    end
end
D(5).y = D(5).y(:,:,1:n0max,:);
D(5).Ey = D(5).Ey(:,1,1:rmax,:);

D(6).name = 'NGF dose response - TrkApErk - ColI';
D(6).n_dim = 2;
D(6).measurand = {'rel. TrkA conc.', 'rel. pErk conc.'};
D(6).t = 60;
D(6).u = [0,0.008,0.04,0.2,1,5];
D(6).y = nan(length(D(6).u),1,10000,2);
D(6).Ey = nan(length(D(6).u),1,10,2);
for r = 1:length(ExpC(1).replicate)
	D(6).replicate{r} = nan(length(D(6).u),1,10000,2);
end

n0max = 0;
rmax = 0;
for d = 1:length(D(6).u)
    n0 = 0;
    for r = 1:4
        nx = size(ExpC(d).replicate(r).data{2},1);
        D(6).y(d,1,n0+1:n0+nx,:) = s3(r)*ExpC(d).replicate(r).data{2}(1:nx,:);
        D(6).Ey(d,1,r,:) = s3(r)*mean(ExpC(d).replicate(r).data{2}(1:nx,:));
        D(6).replicate{r}(d,1,1:nx,:) = s3(r)*ExpC(d).replicate(r).data{2}(1:nx,:);
        n0 = n0 + nx;
    end
        if n0>n0max
      n0max = n0;
    end
    if r > rmax
        rmax = r;
    end
end
D(6).y = D(6).y(:,:,1:n0max,:);
D(6).Ey = D(6).Ey(:,1,1:rmax,:);


ExpC = ExpC4;
D(7).name = 'NGF dose response - ErkpErk - PDL';
D(7).n_dim = 2;
D(7).measurand = {'rel. Erk conc.', 'rel. pErk conc.'};
D(7).t = 60;
D(7).u = [0,0.008,0.04,0.2,1,5];
D(7).y = nan(length(D(7).u),1,10000,2);
D(7).Ey = nan(length(D(7).u),1,10,2);
for r = 1:length(ExpC(1).replicate)
	D(7).replicate{r} = nan(length(D(7).u),1,10000,2);
end


n0max = 0;
rmax = 0;
for d = 1:length(D(7).u)
    n0 = 0;
    for r = 1:4
        nx = size(ExpC(d).replicate(r).data{1},1);
        D(7).y(d,1,n0+1:n0+nx,:) = s4(r)*ExpC(d).replicate(r).data{1}(1:nx,:);
        D(7).Ey(d,1,r,:) = s4(r)*mean(ExpC(d).replicate(r).data{1}(1:nx,:));
        D(7).replicate{r}(d,1,1:nx,:) = s4(r)*ExpC(d).replicate(r).data{1}(1:nx,:);
        n0 = n0 + nx;
    end
    if n0>n0max
      n0max = n0;
    end
    if r > rmax
        rmax = r;
    end
end
D(7).y = D(7).y(:,:,1:n0max,:);
D(7).Ey = D(7).Ey(:,1,1:rmax,:);



D(8).name = 'NGF dose response - ErkpErk - ColI';
D(8).n_dim = 2;
D(8).measurand = {'rel. Erk conc.', 'rel. pErk conc.'};
D(8).t = 60;
D(8).u = [0,0.008,0.04,0.2,1,5];
D(8).y = nan(length(D(8).u),1,10000,2);
D(8).Ey = nan(length(D(8).u),1,10,2);
for r = 1:length(ExpC(1).replicate)
	D(8).replicate{r} = nan(length(D(8).u),1,10000,2);
end

n0max = 0;
rmax = 0;
for d = 1:length(D(8).u)
    n0 = 0;
    for r = 1:4
        nx = size(ExpC(d).replicate(r).data{2},1);
        D(8).y(d,1,n0+1:n0+nx,:) = s4(r)*ExpC(d).replicate(r).data{2}(1:nx,:);
        D(8).Ey(d,1,r,:) = s4(r)*mean(ExpC(d).replicate(r).data{2}(1:nx,:));
        D(8).replicate{r}(d,1,1:nx,:) = s4(r)*ExpC(d).replicate(r).data{2}(1:nx,:);
        n0 = n0 + nx;
    end
  
    if n0>n0max
      n0max = n0;
    end
    if r > rmax
        rmax = r;
    end
end
D(8).y = D(8).y(:,:,1:n0max,:);
D(8).Ey = D(8).Ey(:,1,1:rmax,:);

for e=1:8
    D(e).y(find(D(e).y<0)) = NaN;
end
save data_PDL_ColI D
