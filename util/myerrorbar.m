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