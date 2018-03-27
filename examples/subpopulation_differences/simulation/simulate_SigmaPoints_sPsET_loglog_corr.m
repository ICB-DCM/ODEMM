%function [r1,r2,r3,y,r4,sy] = simulate_SigmaPoints_sPsET_log(t,xi,u)
function varargout = simulate_SigmaPoints_sPsET_loglog_corr(t,xi,u)

%function [y,sy] = simulate_SigmaPoints_sPsET_log(t,xi,u)
%r1 = []; r2 = []; r3 = []; r4 = []; 

% xi = log(k1 k2 k4 k5 k3TrkA0 sErk0) vark3TrkA0 varsErk0

% estruct.beta = @(xi) [xi(1);
%                   xi(2);
%                   xi(3);
%                   xi(4);
%                   xi(5);
%                   xi(6)];


estruct.beta = @(xi) [xi(1);
                  xi(2);
                  xi(3);
                  xi(4);
                  xi(5);
                  xi(6)];
estruct.dbetadxi = @(xi) [1,0,0,0,0,0,0,0,0;
            0,1,0,0,0,0,0,0,0;
            0,0,1,0,0,0,0,0,0;
            0,0,0,1,0,0,0,0,0;
            0,0,0,0,1,0,0,0,0;
            0,0,0,0,0,1,0,0,0];
        
        
% D
estruct.delta = @(xi) [xi(7),xi(8),xi(9)];
estruct.ddeltadxi = @(xi)[0,0,0,0,0,0,1,0,0;
             0,0,0,0,0,0,0,1,0;...
             0,0,0,0,0,0,0,0,1];
A = [eye(6)];
B = [zeros(4,2);
    eye(2)];
estruct.phi = @(beta,b) A*beta + B*b;
estruct.dphidbeta = @(beta,b) A;
estruct.dphidb = @(beta,b) B;
estruct.sigma_noise = @(phi) zeros(1,3);
estruct.dsigma_noisedphi = @(phi) zeros(1,3,6);
if(nargout>=5)
    op_SP.nderiv = 1;
else
    op_SP.nderiv = 0;

end
op_SP.req = [1,1,0,0,0];
op_SP.type_D = 'matrix-logarithm';

[status,SP] = getSigmaPointApp_status_mod(@(phi) simulate_ODEmodel_sPsET_loglog(t,phi,u),xi,estruct,op_SP);

if status >= 0 
    y = SP.my;
    n = size(y,2);
    c=n+1;
    if nargout>=5
        sy = SP.dmydxi;
    end
    for i = 1:n
        for j = i:n
         y(:,c) = SP.Cy(:,i,j);
         if nargout>=5
         sy(:,c,:) = squeeze(SP.dCydxi(:,i,j,:));
         end
         c=c+1;
        end
    end
    varargout{1} = status;
    varargout{2} = [];
    varargout{3} = [];
    varargout{4} = y;
    if nargout>=5
        varargout{5} = [];
        varargout{6} = sy;
    end
else
    varargout{1} = status;
    varargout{2} = [];
    varargout{3} = [];
    varargout{4} = [];
    if nargout>=5
        varargout{5} = [];
        varargout{6} = [];
    end
end
% sy = t x observable x xi
end