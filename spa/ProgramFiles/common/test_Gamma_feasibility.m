function [ cond1_norm, cond2_norm, cond3_norm ] = test_Gamma_feasibility( Gamma )
%TEST_GAMMA_FEASIBILITY Summary of this function goes here
%   Detailed explanation goes here

[K,T] = size(Gamma);

cond1 = Gamma'*ones(K,1) - ones(T,1);
cond2 = min(Gamma,zeros(K,T));
cond3 = max(Gamma,ones(K,T)) - ones(K,T);

cond1_norm = norm(cond1,'fro');
cond2_norm = norm(cond2,'fro');
cond3_norm = norm(cond3,'fro');

end

