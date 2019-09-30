function [ X_rec ] = xrec_signal( Theta, Gamma )

[n,K] = size(Theta);
[K,T] = size(Gamma);

% compute reconstructed data
X_rec = zeros(n,T);
for i=1:n
    for k=1:K
        X_rec(i,:) = X_rec(i,:) + Gamma(k,:)*Theta(i,k);
    end
end

end

