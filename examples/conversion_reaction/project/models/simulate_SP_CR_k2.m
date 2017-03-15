function varargout = simulate_SP_CR_k2(t,xi,u)
          
estruct.beta = @(xi) [xi(1);
                  xi(2);
                  xi(3)];
              
estruct.dbetadxi = @(xi) [1,0,0,0;
            0,1,0,0;
            0,0,1,0];
        
estruct.delta = @(xi) [xi(4)];
estruct.ddeltadxi = @(xi)[0,0,0,1];

A = [eye(3)];

B = [0;1;0];

estruct.phi = @(beta,b) A*beta + B*b;
estruct.dphidbeta = @(beta,b) A;
estruct.dphidb = @(beta,b) B;
estruct.sigma_noise = @(phi) zeros(1,1);
estruct.dsigma_noisedphi = @(phi) zeros(1,1,3);
if(nargout>=5)
    op_SP.nderiv = 1;
else
    op_SP.nderiv = 0;

end
op_SP.req = [1,1,0,0,0];
op_SP.type_D = 'diag-matrix-logarithm';
[status,SP] = getSigmaPointApp_status_mod(@(phi) simulate_CR_log(t,phi,u),xi,estruct,op_SP);

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
end