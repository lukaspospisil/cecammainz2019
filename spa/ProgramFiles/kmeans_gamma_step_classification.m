function [Gamma] = kmeans_gamma_step_classification(X,S)

T = size(X,2);
K = size(S,2);

Gamma = zeros(K,T);
for k=1:K
    Gamma(k,:) = sum((X - S(:,k)).^2,1); % reuse this array
end
[~,idx] = min(Gamma,[],1);
Gamma = zeros(K,T);
Gamma(sub2ind(size(Gamma),idx,1:T)) = 1;

end

