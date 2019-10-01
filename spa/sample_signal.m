clear all

addpath(genpath('ProgramFiles')) % add all files in 'ProgramFiles' directory

%% define sample signal
Trepeat = 1;
n = 3;

[X_orig, Theta_true] = generate_signal_n_K3(Trepeat,n);

% add noise
sigma = 1.5;
X = X_orig + sigma*randn(size(X_orig));

disp(['Signal: [n,T] = ' num2str(size(X,1)) ' x ' num2str(size(X,2)) ' = ' num2str(prod(numel(X)))]);

%% solve clustering problem
% set options of clustering algorithm
in.X = X;
in.options = ClusteringOptions();
in.options.type = 'FEMH1'; % 'FEMH1'/'SPAH1'
in.options.K = 3; % set number of clusters
in.options.epssqr = 1e-1;%1e3; % regularization parameter
in.options.Theta_given = [];%[1,2,3]; % given parameters of clusters
in.options.qp_eps = 1e-3;
in.options.qp_maxit = 2e3;
in.options.eps = 1e-3;%size(X,2)*1e-3;
in.options.maxit = 1e2;
in.options.nanneal = 10; % number of annealing steps

in.options.Sreg = 1e-2;%5e0;%10;%1e-4; % regularization of S problem (only for SPA)
in.options.Sbox = [];%[-5,5]; % constraints S_i \in [min,max], Sbox=[] then no constraints (only for SPA)

in.options.debug = false; % compute and print additional info about progress

% solve the problem using clustering method
tic_clustering = tic;
[out] = signalclustering(in);
time_clustering = toc(tic_clustering);
disp(['Problem solved in ' num2str(time_clustering) 's'])

% compute reconstructed signal
[ X_rec ] = xrec_signal( out.Theta, out.Gamma );

% plot reconstruction
Ttoplot = 1:min(1e4,size(X_rec,2));
ntoplot = 1:size(X,1);
if true
    plot_signal( X_orig(ntoplot,Ttoplot), X(ntoplot,Ttoplot), X_rec(ntoplot,Ttoplot) );
end

% plot Gamma
if true
    plot_signal_clustering( out.Gamma(:,Ttoplot), out.Theta );
end


