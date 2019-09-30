function DG = get_laplacian_grid3d(T,R2,R1)

DR = get_laplacian_grid2d(R2,R1); % space regularization
DT = get_laplacian_grid1d( T ); % time regularization

DG = kron(speye(T),DR) + kron(DT,speye(R1*R2));

end
