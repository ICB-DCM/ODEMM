function s = getScalingFactors(varargin)
% Calculates scaling factors for replicates such that the distance between
% the means in log-space are minimal.
%
% USAGE:
% s = getScalingFactors(ExpC)\n
% s = getScalingFactors(ExpC_1,ExpC_2)
%
% Parameters:
% varargin:
% ExpC: struct of experiments
%
% Required fields of ExpC:
%   name: string specifying the conditions
%   time: time point of measurement
%   stimulus: stimulus for measurement
%   replicate: struct of replicates
%   * name: string specifying the replicate
%   * measurands: names of measurands
%   * ndata: matrices under different conditions
%         (one row represents one observed cell with the
%          data in the order of the measurands. The different
%          rows provide measurement data for different cells)
%
% Return values:
% s: (1 x n_r) vector including scaling factor for every replicate


%%  extract the different time and stimulus conditions (u_t) and the total number of replicate n_rtot
n_rtot = 0;
n_exp = nargin;
u_t = [];
for e = 1:n_exp
    ExpC = varargin{e};
    n_rtot = n_rtot + length(ExpC(1).replicate);
    for j = 1:length(ExpC)
        if isempty(u_t) || ~ismember([ExpC(j).stimulus,ExpC(j).time],u_t,'rows')
            u_t = [u_t; [ExpC(j).stimulus,ExpC(j).time]];
        end
    end
end
n_c = length(ExpC(1).replicate(1).data);
n_diff = size(u_t,1)*n_c;
n_meas = length(ExpC(1).replicate(1).measurands);
M = nan(n_diff,n_meas,n_rtot);
r_count = 0;
for e = 1:n_exp
    ExpC = varargin{e};
    for j = 1:length(ExpC)
        [~,c] = ismember([ExpC(j).stimulus,ExpC(j).time],u_t,'rows');
        for r = 1:length(ExpC(j).replicate)
            if ~iscell(ExpC(j).replicate(r).data)
                M(c,r_count+r) = mean(ExpC(j).replicate(r).data(:,1)');
            else
                for b = 1:length(ExpC(j).replicate(r).data)
                    nx = size(ExpC(j).replicate(r).data{b},1);
                    M(c+(size(u_t,1)*(b-1)),:,r_count+r) = mean(ExpC(j).replicate(r).data{b}(1:nx,:));
                end
            end
        end
    end
    r_count = length(ExpC(j).replicate);
end

scale = 'log';

switch scale
    case 'log'
        switch n_meas
            case 2
                J2 = @(s) nansum(nansum(bsxfun(@minus,(bsxfun(@plus,...
                    log(squeeze(M(:,1,:))),log(s)')),...
                    (nanmean(bsxfun(@plus,log(squeeze(M(:,1,:))),log(s)'),2))).^2))+...
                    nansum(nansum(bsxfun(@minus,(bsxfun(@plus,...
                    log(squeeze(M(:,2,:))),log(s)')),...
                    (nanmean(bsxfun(@plus,log(squeeze(M(:,2,:))),log(s)'),2))).^2));
            case 1
                J2 = @(s) nansum(nansum(bsxfun(@minus,(bsxfun(@plus,...
                    log(squeeze(M(:,1,:))),log(s)')),...
                    (nanmean(bsxfun(@plus,log(squeeze(M(:,1,:))),log(s)'),2))).^2));
        end
    case 'lin'
        switch n_meas
            case 2
                J2 = @(s) nansum(nansum(bsxfun(@minus,(bsxfun(@times,...
                    (squeeze(M(:,1,:))),(s)')),...
                    (nanmean(bsxfun(@times,(squeeze(M(:,1,:))),(s)'),2))).^2))+...
                    nansum(nansum(bsxfun(@minus,(bsxfun(@times,...
                    (squeeze(M(:,2,:))),(s)')),...
                    (nanmean(bsxfun(@times,(squeeze(M(:,2,:))),(s)'),2))).^2));
            case 1
                J2 = @(s) nansum(nansum(bsxfun(@minus,(bsxfun(@times,...
                    (squeeze(M(:,1,:))),(s)')),...
                    (nanmean(bsxfun(@times,(squeeze(M(:,1,:))),(s)'),2))).^2));
        end
end

s = fmincon(@(s) J2(s),ones(n_rtot,1),[],[],ones(1,n_rtot),n_rtot);

