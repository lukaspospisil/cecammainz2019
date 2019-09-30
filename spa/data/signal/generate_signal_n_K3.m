function [ X_true, Theta_true ] = generate_signal_n_K3( Trepeat, n )

% affiliation for one period
aff_one_period = [1*ones(1,100) 2*ones(1,80) 1*ones(1,50) 3*ones(1,70) 2*ones(1,90) 1*ones(1,50) 3*ones(1,60)];
K = 3;

T = length(aff_one_period);
Gamma = zeros(K,T);
for k=1:K
    Gamma(k,aff_one_period == k) = 1;
end

Theta_true = [1, 2, -1, 0; ...
         2, -1, 0, 1; ...
         -1, 2, 0, 1.5]';

if n > 4
    Theta_true = [Theta_true;zeros(n-4,K)];
end

X_true = kron(ones(1,Trepeat),Theta_true*Gamma);

end

