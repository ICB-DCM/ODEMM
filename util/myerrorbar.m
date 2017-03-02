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

function  [] = myerrorbar(X,Y,E,xlength)
    errorbar(X,Y,E,'Color','k');
    for k = 1:length(X)
     x = [X(k) - xlength, X(k) + xlength];
     y_h = [Y(k) + E(k), Y(k) + E(k)];
     line(x, y_h,'Color','k');
     y_b = [Y(k) - E(k), Y(k) - E(k)];
     line(x, y_b,'Color','k');
    end
end