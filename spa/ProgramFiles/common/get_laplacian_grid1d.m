function [ A ] = get_laplacian_grid1d( n )

if n > 1
    ondiag = [ 1 2*ones(1,n-2) 1];
    offdiag = -ones(1,n);
    A = spdiags([offdiag' ondiag' offdiag'],-1:1,n,n);
else
    A = 0;
end

end

