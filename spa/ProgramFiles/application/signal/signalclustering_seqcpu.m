function [ out ] = signalclustering_seqcpu( in )

% in.X in format [n,T]
[n,T] = size(in.X);

Hreg = get_laplacian_grid1d(T); % time regularization matrix

tic_solver = tic;
if strcmp(in.options.type,'FEMH1')
    [ out.Theta, out.Gamma, out.L, out.Llin, out.Lquad ] = femH1(in.X, Hreg, in.options);
end
if strcmp(in.options.type,'SPAH1')
    [ out.Theta, out.Gamma, out.L, out.Llin, out.Lquad ] = spaH1(in.X, Hreg, in.options);
end

time_solver = toc(tic_solver);
if in.options.debug
    disp(['Pure solver time: ' num2str(time_solver) 's'])
end

end