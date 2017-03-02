function [] = test_onestage_SP(t,xi,u)
                
estruct.beta = @(xi) [xi(1);
                  xi(2);
                  xi(3);
                  xi(4);
                  xi(4)];
              
estruct.dbetadxi = @(xi) [1,0,0,0,0;
                          0,1,0,0,0;
                          0,0,1,0,0;
                          0,0,0,1,0;
                          0,0,0,1,0];
        
        
% D
estruct.delta = @(xi) [xi(5),xi(5)];
estruct.ddeltadxi = @(xi)[0,0,0,0,1;
                         0,0,0,0,1];

A = [eye(5)];

B = [0,0;
     0,0;
     0,0;
     1,0;
     0,1];

estruct.phi = @(beta,b) A*beta + B*b;
estruct.dphidbeta = @(beta,b) A;
estruct.dphidb = @(beta,b) B;
estruct.sigma_noise = @(phi) zeros(1,1);
estruct.dsigma_noisedphi = @(phi) zeros(1,1,1);
if(nargout>=5)
    op_SP.nderiv = 1;
else
    op_SP.nderiv = 0;

end
op_SP.req = [1,1,0,0,0];
op_SP.type_D = 'diag-matrix-logarithm';
[SP] = testSigmaPointApp_mod(@(phi) simulate_oneStage_RRE(t,phi,u),xi,estruct,op_SP);

end