% Plotting routine for errorbars. Plots X versus Y with errorbars of size E
% and length xlength.
%
% USAGE:
% myerrorbar(X,Y,E,xlength)
%
% Parameters: 
% X: Values at x-axis
% Y: Values at y-axis
% E: symmetric errorbars of size 2*E
% xlength: width of errorbar

function  [] = myerrorbar(varargin)
if nargin <= 4
    col = 'k';
else
    col = varargin{5};
end

X = varargin{1};
Y = varargin{2};
E = varargin{3};
xlength = varargin{4};

if nargin < 6
    marker = '.';
    markersize = 8;
else
    marker = varargin{6};
    markersize = varargin{7};
end
    for k = 1:length(X)
     x = [X(k) - xlength, X(k) + xlength];
     y_h = [Y(k) + E(k), Y(k) + E(k)];
     line(x, y_h,'Color',col);hold on;
     y_b = [Y(k) - E(k), Y(k) - E(k)];
     line(x, y_b,'Color',col);
     line([X(k),X(k)],[Y(k)+E(k),Y(k)-E(k)],'Color',col);
     if nargin > 5
        plot(X(k),Y(k),marker,'MarkerSize',markersize,'color',col);
     end
    end
end