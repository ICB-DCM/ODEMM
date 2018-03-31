function [conditions,D] = collectConditions(D,M)
% This function collects all different conditions regarding input/differences
%   between subpopulations/experiments and timepoints.
%
% USAGE:
% [conditions,D] = collectConditions(D,M)
%
% Parameters:
%   D: data struct
%   M: model struct
%
% Return values:
%   conditions: conditions struct
%   D: updated data struct
%
% Required fields of D:
%    t: time vector
%    u: vector of stimulations
%
% Required fields of M:
%    n_subpop: number of subpopulations
%    u{s,e}: input vector capturing differences between subpopulations
%               and experiments
%
% Generated fields of D:
%     c: n_subpop x (n_u + n_differences) matrix linking condition to
%           data
%
% Generated fields of conditions: 
%      input: (n_u + n_differences) x 1 input vector
%      time: 1 x n_t time vector
%      sigma: 1 x n_t vector of sigmas for condition c


conditions(1).input = [D(1).u(:,1);M.u{1,1}];
conditions(1).time = D(1).t;
D(1).c = nan(M.n_subpop,size(D(1).u,2));
D(1).c(1,1) = 1;
for e = 1:length(D)
    for s=1:M.n_subpop
        for d = 1:size(D(e).u,2)
            for k = 1:size(D(e).t,1)
                add = true;
                for c = 1:length(conditions)
                    if ~isequal(conditions(c).input,[D(e).u(:,d);M.u{s,e}])
                        continue;
                    elseif ~ismember(D(e).t(k),conditions(c).time) % dose and theta already there but new time point
                        conditions(c).time = [conditions(c).time, D(e).t(k)];
                        add = false;
                        D(e).c(s,d) = c;
                        break;
                    elseif ismember(D(e).t(k),conditions(c).time)
                        add = false;
                        D(e).c(s,d) = c;
                        break;
                    end
                end
            end
            if add %condition not considered yet -> add
                c = length(conditions)+1;
                conditions(c).input = [D(e).u(:,d);M.u{s,e}];
                conditions(c).time = D(e).t;
                D(e).c(s,d) = c;
            end
        end
    end
end

% ensure increasing time vector
for c = 1:length(conditions)
    conditions(c).time = sort(conditions(c).time);
end
end

