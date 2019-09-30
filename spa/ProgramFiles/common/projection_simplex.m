function PX = projection_simplex(X,K)

T = length(X)/K;

XX = reshape(X,T,K)';

PXX = projection_simplexes(XX);

PX = reshape(PXX',K*T,1);

end