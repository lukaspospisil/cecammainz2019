function [Gamma] = spa_gamma_step_classification(X,S)

[n,T] = size(X);
K = size(S,2);

% compute new gamma
HS=S'*S; HS=(HS+HS');
H = kron(HS,speye(T));
CTX = 2*reshape((X'*S),K*T,1);
        
normH = gershgorin(H);
        
% solve qp
[gamma_vec] = spgqp(H, CTX, rand(size(CTX)), K, normH, 1e-6, 1e3);
        
Gamma = reshape(gamma_vec,T,K)';

end

