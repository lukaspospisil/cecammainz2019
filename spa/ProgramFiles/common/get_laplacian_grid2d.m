function DG = get_laplacian_grid2d(R2,R1)

DG = (kron(speye(R2),get_laplacian_grid1d(R1)) + kron(get_laplacian_grid1d(R2),speye(R1)));

end
