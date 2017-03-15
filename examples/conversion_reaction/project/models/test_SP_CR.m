function [] = test_SP_CR(t,xi,u)
% function to test the Sigma Point Approximation for the conversion
% reaction for the ground truth model.
%
% Parameters:
% t: time vector
% xi: parameter vector
% u: simulation
             
estruct.beta = @(xi) [xi(1);
                  xi(2);
                  xi(3)];
              
estruct.dbetadxi = @(xi) [1,0,0,0;
            0,1,0,0;
            0,0,1,0];
        
estruct.delta = @(xi) [xi(4)];
estruct.ddeltadxi = @(xi)[0,0,0,1];

A = [eye(3)];

B = [0;0;1];

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

[SP] = testSigmaPointApp(@(phi) simulate_CR_log(t,phi,u),xi,estruct,op_SP);
end