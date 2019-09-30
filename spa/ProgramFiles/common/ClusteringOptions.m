classdef ClusteringOptions
    %CLUSTERING_OPTIONS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nanneal         % number of annealing steps
        
        type            % 'FEMH1'/'SPAH1'/'DDFEMH1'/'DDSPAH1'
        K               % number of clusters
        epssqr          % regularization parameter

        Sreg            % regularization of S problem (only for SPA)
        Sbox            % constraints S_i \in [min,max], Sbox=[] then no constraints (only for SPA)

        m               % parameter of SPAS method
        
        maxit           % max number of outer (subspace) algorithm (only if Theta is not given)
        eps             % stopping criteria of outer (subspace) algorithm (only if Theta is not given)
        qp_maxit        % max number of QP solver iterations
        qp_eps          % stopping criteria of QP solver
        
        Theta_given     % given parameters of clusters
        Gamma0          % given initial approximation of subspace algorithm
        
        debug           % compute and display some additional informations
    end
    
    methods
        function obj = ClusteringOptions()
            % constructor - set default values
            obj.nanneal = 1;
            
            obj.type = 'FEMH1';
            obj.K = 2;
            obj.epssqr = 1e-2;
            
            obj.Sreg = 1e-5;
            obj.Sbox = [];
            
            obj.m = 10;
            
            obj.maxit = 100;
            obj.eps = 1e-4;
            obj.qp_maxit = 1e3;
            obj.qp_eps = 1e-4;
            
            obj.Theta_given = [];
            obj.Gamma0 = [];
            
            obj.debug = false;
        end
        
    end
end    
