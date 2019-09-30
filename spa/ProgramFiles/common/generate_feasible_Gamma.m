function [ Gamma ] = generate_feasible_Gamma( K,T,Gamma0 )

% generate initial random feasible Gamma [K,T]
if ~isempty(Gamma0)
    Gamma = Gamma0;
else
    Gamma = rand(K,T);
end

% initial gamma has to be feasible (i.e sum to one for each t,r and positive)
Gamma_sum = sum(Gamma,1);
for k=1:K
    Gamma(k,:) = Gamma(k,:)./Gamma_sum;
end

end

