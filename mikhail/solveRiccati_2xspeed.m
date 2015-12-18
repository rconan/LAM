function [gamma,tr,condnumber] = solveRiccati_2Xspeed(A,C,Q,R, MaxIter);
% ------------HEADER-----------------
% Objective         ::  Solve the Riccati equation using a doubling speed
%                       algorithm - See Optimal Filtering from
%                       Anderson&Moore, pag. 159
% INPUT VARS
% A                 ::  State transition matrix (Nstates x Nstates)
% C                 ::  Measurement matrix (Nmeasurements x Nstates)
% Q                 ::  State error covariance matrix (Nstates x Nstates)
% R                 ::  Measurement covariance matrix (Nmeasurements x Nmeasurements)
% MaxIter           ::  Number of maximum iterations
%
% OUTPUT VARS
% gamma             ::  The solution of the Riccati equation of order
%                       Sigma_{2^iter | 2^iter-1}
%
%Created by         ::  Carlos Correia
%Creation date      ::  August 08
%Change Record:     ::  
% ------------HEADER END----------------


% --- parse input ---
if nargin < 5
    MaxIter = 50;
end

% --- initialisation ---
alpha   = A';
beta    = C'*inv(R)*C;
gamma   = Q;
Id = eye(size(Q));
for iter = 1:MaxIter
    common  = inv(Id + beta*gamma);
    alpha1   = alpha*common*alpha;
    beta1    = beta + alpha*common*beta*alpha';
    gamma1   = gamma + alpha'*gamma*common*alpha;
    
    if nargout > 1
        tr(iter) = trace(gamma1 - gamma);
    end
    if nargout > 2
        condnumber(iter) = cond(gamma1);
    end
    
    alpha = alpha1;
    beta = beta1;
    gamma = gamma1;
    
end

    
    


