function [] = printParams(parameters,varargin)
% Help function to print parameters names and values.
%
% USAGE:
% [] = printParams(parameters,xi)
%
% Parameters:
% parameters: parameters struct
% varargin:
% xi: parameter values printed along with parameter name
%
% Required fields of parameters:
% name = {'name1',...} names of parameters
%
if nargin >= 2
    xi = varargin{1};
end
if nargin >=2
    for i = 1:parameters.number
        disp([num2str(i) ' ' parameters.name{i} '   ' num2str(xi(i))]);
    end
else
    for i = 1:parameters.number
        disp([num2str(i) ' ' parameters.name{i}]);
    end
end