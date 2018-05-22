function [parameters, conditions, varargout] = getRRErhos(parameters,conditions, varargin)
% This function defines the parameters needed for the parametrization of
% the variances in case Reaction Rate Equations are used for the
% mechanistic description of the means
%
% USAGE:
% [parameters,conditions] = getRRErhos(parameters,conditions) \n
% [parameters,conditions,D] = getRRErhos(parameters,conditions,options,D,M)
%
% Parameters:
% parameters: parameters struct
% conditions: conditions struct (see collectConditions.m)
% varargin:
%   options:
%   D: data struct (see logLikelihood.m)
%   M: model struct (see generateODEMM.m)
%
% Return values:
% parameters: updated parameters struct
% conditions: updated conditions struct
% D: updated data struct
%
%  Optional fields of options:
%  rhos: \n
%      = ''time-dependent'': (default) one rho for every subpopulation and time point \n
%      = ''only-one'': only one rho for everything \n
%      = ''subpopulation-specific'': for every subpopulation one rho \n
%  boundaries: boundaries for optimization for the rho parameters with
%  fields
%    * min
%    * max
%
% Generated fields of parameters:
%   names: names for rho parameters are added
%
% Generated fields of D:
%   rho: (if n_dim = 1)
%   rho: (if n_dim = 2)
%

options.rhos = 'time-dependent';
if nargin >= 3
    options = varargin{1};
end
if nargin >=4
    D = varargin{2};
end
if nargin >=5
    M = varargin{3};
end
disp('Start assigning rho parameters for every conditions...')
num_params = length(parameters.name);
k = length(parameters.name)+1;

if strcmp(options.rhos,'time-dependent')
    for e = 1:length(D)
        for s = 1:M.n_subpop
            if D(e).n_dim == 1
                D(e).rho{s} = [];
            else
                D(e).rho{s} = [];
            end
            for n_t = 1:length(D(e).t)
                for n_u = 1:size(D(e).u,2)
                    if D(e).n_dim == 1
                        parameters.name{k} = ['\rho_{exp. ' num2str(e) ', t=' num2str(D(e).t(n_t)) ', u=' num2str(D(e).u(1,n_u)), ' s = ' num2str(s) ' }'];
                        D(e).rho{s} = [D(e).rho{s}, k];
                        k = k+1;
                    else
                        if options.covariance
                            D(e).rho{s} = [D(e).rho{s}, [k;k+1;k+2]];
                            parameters.name{k} = ['\rho_{exp}. ' num2str(e) ', t=' num2str(D(e).t(n_t)) ' u=' num2str(D(e).u(1,n_u)) ' s = ' num2str(s) ', var meas.1)}'];
                            parameters.name{k+1} = ['\rho_{exp}. ' num2str(e) ', t=' num2str(D(e).t(n_t)) ' u=' num2str(D(e).u(1,n_u)) ' s = ' num2str(s) ', var meas.2)}'];
                            parameters.name{k+2} = ['\rho_{exp}. ' num2str(e) ', t=' num2str(D(e).t(n_t)) ' u=' num2str(D(e).u(1,n_u)) ' s = ' num2str(s) ', cov)}'];
                            k = k+3;
                        else
                            D(e).rho{s} = [D(e).rho{s}, [k;k+1]];
                            parameters.name{k} = ['\rho_{exp}. ' num2str(e) ', t=' num2str(D(e).t(n_t)) ' u=' num2str(D(e).u(1,n_u)) ' s = ' num2str(s) ', var meas.1)}'];
                            parameters.name{k+1} = ['\rho_{exp}. ' num2str(e) ', t=' num2str(D(e).t(n_t)) ' u=' num2str(D(e).u(1,n_u)) ' s = ' num2str(s) ', var meas.2)}'];
                            k = k+2;
                        end
                    end
                end
            end
        end
    end
    
elseif strcmp(options.rhos,'subpopulation-specific')
    
    for e = 1:length(D)
        for s = 1:2
            if D(e).n_dim == 1
                D(e).rho{s} = [];
            else
                D(e).rho{s} = [];
            end
            if D(e).n_dim == 1
                parameters.name{k} = ['\rho_{exp. ' num2str(e) ',  s = ' num2str(s) ' }'];
                D(e).rho{s} = [D(e).rho{s}, k];
                k = k+1;
            else
                if options.covariance
                    D(e).rho{s} = [D(e).rho{s}, [k;k+1;k+2]];
                    parameters.name{k} = ['\rho_{exp}. ' num2str(e) ' s = ' num2str(s) ', var meas.1)}'];
                    parameters.name{k+1} = ['\rho_{exp}. ' num2str(e) ' s = ' num2str(s) ', var meas.2)}'];
                    parameters.name{k+2} = ['\rho_{exp}. ' num2str(e) ' s = ' num2str(s) ', cov)}'];
                    k = k+3;
                else
                    D(e).rho{s} = [D(e).rho{s}, [k;k+1]];
                    parameters.name{k} = ['\rho_{exp}. ' num2str(e) ' s = ' num2str(s) ', var meas.1)}'];
                    parameters.name{k+1} = ['\rho_{exp}. ' num2str(e) ' s = ' num2str(s) ', var meas.2)}'];
                    k = k+2;
                end
            end
        end
    end
else
    for c=1:length(conditions)
        switch options.dimension
            case 'univariate'
                conditions(c).rho = [];
            case 'multivariate'
                conditions(c).rho = [];
        end
        switch options.rhos
%             case 'condition-dependent'
%                 for n_t = 1:length(conditions(c).time)
%                     parameters.name{k} = ['\rho_{cond. ' num2str(c) ', t=' num2str(conditions(c).time(n_t)) ')}'];
%                     conditions(c).rho = [conditions(c).rho, k];
%                     k = k+1;
%                 end
            case 'time-independent'
                parameters.name{k} = ['\rho_{cond. ' num2str(c) ')}'];
                conditions(c).rho = [conditions(c).rho, k];
                k = k+1;
            case 'only-one'
                if c == 1
                    parameters.name{k} = ['\rho'];
                end
                conditions(c).rho = [k];
                if c == numel(conditions)
                    k=k+1;
                end
        end
    end
end
parameters.number = length(parameters.name);
if ~isfield(options,'boundaries')
    options.boundaries.min = 0;
    options.boundaries.max = 1;
end
parameters.max = [parameters.max; options.boundaries.max*ones(k-num_params-1,1)];
parameters.min = [parameters.min; options.boundaries.min*ones(k-num_params-1,1)];

disp('...done!')
disp([num2str(k-num_params-1) ' rhos have been introduced'])
if nargout >= 1
    varargout{1} = D;
end
