clear all

addpath(genpath('ProgramFiles')) % add all files in 'ProgramFiles' directory

T = 1e2; % the length of time-series

% exact discretization parameters (mean value)
SX_exact = [1,2,3];
SY_exact = [4,5];

KX = length(SX_exact);
KY = length(SY_exact);

Lambda_exact = [1, 0, 1; 0 1 0];

% generate random signal X
Xidx = randi(length(SX_exact),1,T);
X_exact = zeros(1,T);
GammaX_exact = zeros(length(SX_exact),T);

for i=1:KX
    X_exact(Xidx == i) = SX_exact(i);
    GammaX_exact(i,Xidx == i) = 1;
end

GammaY_exact = Lambda_exact*GammaX_exact;
Y_exact = SY_exact*GammaY_exact;

if false
    figure
    subplot(2,1,1)
    hold on
    plot(1:T,X_exact,'blue');
    xlabel('$t$','interpreter','latex')
    ylabel('$X_t$','interpreter','latex')
    hold off
    
    subplot(2,1,2)
    hold on
    plot(1:T,Y_exact,'red');
    xlabel('$t$','interpreter','latex')
    ylabel('$Y_t$','interpreter','latex')
    hold off
end

%
% TEST 1:
% - take X_exact, Y_exact
% - perform the classification using SX_exact, SY_exact ( compute GammaX, GammaY using K-means)
%   (copy-paste Gamma step from K-means)
% - reconstruct Lambda from "classic" markov GammaY = Lambda*GammaX
%   (see KL + Jensen estimation)
%
GammaX = kmeans_gamma_step_classification(X_exact,SX_exact);
GammaY = kmeans_gamma_step_classification(Y_exact,SY_exact);
Lambda1 = markov_est(GammaX,GammaY);

% add noise to the signal
sigma = 0.3;
X = X_exact + sigma*rand(size(X_exact));
Y = Y_exact + sigma*rand(size(Y_exact));

if false
    figure
    subplot(2,1,1)
    hold on
    plot(1:T,X,'blue');
    xlabel('$t$','interpreter','latex')
    ylabel('$X_t$','interpreter','latex')
    hold off
    
    subplot(2,1,2)
    hold on
    plot(1:T,Y,'red');
    xlabel('$t$','interpreter','latex')
    ylabel('$Y_t$','interpreter','latex')
    hold off
end

%
% TEST 2:
% - throw away X_exact, Y_exact and work only with "continuous" X,Y
% - perform the classification using SX_exact, SY_exact ( compute GammaX, GammaY using K-means)
%   (copy-paste Gamma step from K-means) = binary Gamma
% - reconstruct Lambda from "classic" markov GammaY = Lambda*GammaX
%   (see KL + Jensen estimation)
%
GammaX = kmeans_gamma_step_classification(X,SX_exact);
GammaY = kmeans_gamma_step_classification(Y,SY_exact);
Lambda2 = markov_est(GammaX,GammaY);

markov_eps = 1;

%% solve clustering problem
% set options of clustering algorithm
in.X = [Y;markov_eps*X];
in.options = ClusteringOptions();
in.options.type = 'SPAH1'; % 'FEMH1'/'SPAH1'
in.options.K = 3; % set number of clusters
in.options.epssqr = 0;%1e3; % regularization parameter
in.options.Theta_given = [];%[1,2,3]; % given parameters of clusters
in.options.qp_eps = 1e-6;
in.options.qp_maxit = 5e3;
in.options.eps = 1e-4;%size(X,2)*1e-3;
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

% process results back to original setting
SX = out.Theta(2,:)/markov_eps;
SY_Lambda = out.Theta(1,:);
GammaX = out.Gamma;


