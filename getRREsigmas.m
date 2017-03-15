function [parameters, conditions, varargout] = getRREsigmas(parameters,conditions, varargin)
% This function defines the parameters needed for the parametrization of
% the variances in case Reaction Rate Equations are used for the
% mechanistic description of the means
%
% USAGE:
% [parameters,conditions] = getRREsigmas(parameters,conditions) \n
% [parameters,conditions,D] = getRREsigmas(parameters,conditions,options,D,M)
%
% Parameters:
% parameters: parameters struct 
% conditions: conditions struct (see ...)
% varargin:
%   options:
%   D: data struct (see logLikelihood.m)
%   M: model struct (see generateODEMM.m)
%
% Return values:
% parameters: updated parameters struct
% conditions: updated conditions struct
% varargout:
%   D: updated data struct
%
%
%  Optional fields of options:
%  sigmas: \n
%         = 'condition-dependent': (default) assign sigma for every time point \n
%         = 'time-dependent': one sigma for every subpopulation and time point \n
%         = 'only-one': only one sigma for everything \n
%         = 'subpopulation-specific': for every subpopulation one sigma \n
%  boundaries
%    * .min
%    * .max
% 
% Generated fields of parameters:
%   names: names for sigma parameters are added
%   
% Generated fields of D:
%   sigma: (if n_dim = 1)
%   Sigma: (if n_dim = 2)
%   
%

options.sigmas = 'condition-dependent';
if nargin >= 3
    options = varargin{1};
end
if nargin >=4 
    D = varargin{2};
end
if nargin >=5
    M = varargin{3};
end
disp('Start assigning sigma parameters for every conditions...')
num_params = length(parameters.name);
k = length(parameters.name)+1;

if strcmp(options.sigmas,'time-dependent') 
    for e = 1:length(D)
        for s = 1:M.n_subpop
             if D(e).n_dim == 1
                D(e).sigma{s} = [];
             else
                D(e).Sigma{s} = [];
            end
            for n_t = 1:length(D(e).t)
               for n_u = 1:size(D(e).u,2)
               if D(e).n_dim == 1
                    parameters.name{k} = ['log_{10}(\sigma_{exp. ' num2str(e) ', t=' num2str(D(e).t(n_t)) ', u=' num2str(D(e).u(1,n_u)), ' s = ' num2str(s) ' }'];
                    D(e).sigma{s} = [D(e).sigma{s}, k];
                    k = k+1;        
               else
                   if options.covariance
                        D(e).Sigma{s} = [D(e).Sigma{s}, [k;k+1;k+2]];
                        parameters.name{k} = ['log_{10}(\Sigma_{exp}. ' num2str(e) ', t=' num2str(D(e).t(n_t)) ' u=' num2str(D(e).u(1,n_u)) ' s = ' num2str(s) ', var meas.1)}'];
                        parameters.name{k+1} = ['log_{10}(\Sigma_{exp}. ' num2str(e) ', t=' num2str(D(e).t(n_t)) ' u=' num2str(D(e).u(1,n_u)) ' s = ' num2str(s) ', var meas.2)}'];
                        parameters.name{k+2} = ['log_{10}(\Sigma_{exp}. ' num2str(e) ', t=' num2str(D(e).t(n_t)) ' u=' num2str(D(e).u(1,n_u)) ' s = ' num2str(s) ', cov)}'];
                        k = k+3;    
                   else
                        D(e).Sigma{s} = [D(e).Sigma{s}, [k;k+1]];
                        parameters.name{k} = ['log_{10}(\Sigma_{exp}. ' num2str(e) ', t=' num2str(D(e).t(n_t)) ' u=' num2str(D(e).u(1,n_u)) ' s = ' num2str(s) ', var meas.1)}'];
                        parameters.name{k+1} = ['log_{10}(\Sigma_{exp}. ' num2str(e) ', t=' num2str(D(e).t(n_t)) ' u=' num2str(D(e).u(1,n_u)) ' s = ' num2str(s) ', var meas.2)}'];
                        k = k+2;   
                   end
                end
               end
            end
        end
    end
    
elseif strcmp(options.sigmas,'subpopulation-specific')
  
    for e = 1:length(D)
        for s = 1:2
             if D(e).n_dim == 1
                D(e).sigma{s} = [];
             else
                D(e).Sigma{s} = [];
             end
               if D(e).n_dim == 1
                    parameters.name{k} = ['log_{10}(\sigma_{exp. ' num2str(e) ',  s = ' num2str(s) ' }'];
                    D(e).sigma{s} = [D(e).sigma{s}, k];
                    k = k+1;        
               else
                   if options.covariance
                        D(e).Sigma{s} = [D(e).Sigma{s}, [k;k+1;k+2]];
                        parameters.name{k} = ['log_{10}(\Sigma_{exp}. ' num2str(e) ' s = ' num2str(s) ', var meas.1)}'];
                        parameters.name{k+1} = ['log_{10}(\Sigma_{exp}. ' num2str(e) ' s = ' num2str(s) ', var meas.2)}'];
                        parameters.name{k+2} = ['log_{10}(\Sigma_{exp}. ' num2str(e) ' s = ' num2str(s) ', cov)}'];
                        k = k+3;    
                   else
                        D(e).Sigma{s} = [D(e).Sigma{s}, [k;k+1]];
                        parameters.name{k} = ['log_{10}(\Sigma_{exp}. ' num2str(e) ' s = ' num2str(s) ', var meas.1)}'];
                        parameters.name{k+1} = ['log_{10}(\Sigma_{exp}. ' num2str(e) ' s = ' num2str(s) ', var meas.2)}'];
                        k = k+2;   
                   end
               end
        end
    end
else
for c=1:length(conditions)
    switch options.dimension
        case 'univariate'
            conditions(c).sigma = [];
        case 'multivariate'
            conditions(c).Sigma = [];
    end
    switch options.sigmas
        case 'condition-dependent'
            for n_t = 1:length(conditions(c).time)
                switch options.dimension
                    case 'univariate'
                    parameters.name{k} = ['log_{10}(\sigma_{cond. ' num2str(c) ', t=' num2str(conditions(c).time(n_t)) ')}'];
                    conditions(c).sigma = [conditions(c).sigma, k];
                    k = k+1;        
                    case 'multivariate'
                    conditions(c).Sigma = [conditions(c).Sigma, [k;k+1;k+2]];    
                    parameters.name{k} = ['log_{10}(\Sigma_{cond. ' num2str(c) ', t=' num2str(conditions(c).time(n_t)) ', var meas.1)}'];
                    parameters.name{k+1} = ['log_{10}(\Sigma_{cond. ' num2str(c) ', t=' num2str(conditions(c).time(n_t)) ', var meas.2)}'];
                    parameters.name{k+2} = ['log_{10}(\Sigma_{cond. ' num2str(c) ', t=' num2str(conditions(c).time(n_t)) ', cov)}'];
                    k = k+3;         
                end
            end
        case 'time-independent'
            parameters.name{k} = ['log_{10}(\sigma_{cond. ' num2str(c) ')}'];
            conditions(c).sigma = [conditions(c).sigma, k];
            k = k+1;   
        case 'only-one'
            if c == 1
               parameters.name{k} = ['log_{10}(\sigma)'];
            end
            conditions(c).sigma = [k];
            if c == numel(conditions)
                k=k+1;
            end
    end  
end
end
parameters.number = length(parameters.name);
if ~isfield(options,'boundaries')
    options.boundaries.min = -2;
    options.boundaries.max = 2;
end
parameters.max = [parameters.max; options.boundaries.max*ones(k-num_params-1,1)];
parameters.min = [parameters.min; options.boundaries.min*ones(k-num_params-1,1)];

disp('...done!')
disp([num2str(k-num_params-1) ' sigmas have been introduced'])
if nargout >= 1
    varargout{1} = D;
end
