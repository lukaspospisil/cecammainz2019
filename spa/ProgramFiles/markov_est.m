function [Lambda] = markov_est(GammaX,GammaY)
% estimate Lambda using KL + jensen (= likelihood)

KX = size(GammaX,1);
KY = size(GammaY,1);

Lambda = zeros(KY,KX);
for kx = 1:KX
    for ky = 1:KY
        Lambda(ky,kx) = dot(GammaX(kx,:),GammaY(ky,:));
    end
    % normalize Lambda
    Lambda(:,kx) = Lambda(:,kx)./sum(Lambda(:,kx));
end

end

