function [ Theta, Gamma, L, Llin, Lquad, Lit ] = femH1(X,Hreg,options )

% provide the signal in format [n,T]
[n,T] = size(X);
K = options.K;

% check sizes of input variables
if size(Hreg,1) ~= T
    error('! size of regularization matrix does not match with the size of data')
end
if ~isempty(options.Theta_given)
    if size(options.Theta_given,1) ~= n
        error('! n of data and size of Theta_given do not match')
    end
    if size(options.Theta_given,2) ~= K
        error('! given K and size of Theta_given do not match')
    end
end
if ~isempty(options.Gamma0)
    if size(options.Gamma0,1) ~= K
        error('! K and size of Gamma0 do not match')
    end
    if size(options.Gamma0,2) ~= T
        error('! length of data and Gamma0 do not match')
    end
end


% if someone provided Theta, then there are not outer iterations
if isempty(options.Theta_given)
    maxit = options.maxit;
else
    maxit = 1;
end

% prepare QP objects
HG = kron(speye(K),Hreg);
normH = gershgorin(Hreg);

g = zeros(K*T,1);

% here will be stored solution of model parameters - one for each cluster
Theta = zeros(n,K);

% initial object function value
L_best = Inf;

for anneal_idx = 1:options.nanneal
    % display some informations, so we know that something is happening
    if options.debug
        disp(['annealing step ' num2str(anneal_idx) ' of ' num2str(options.nanneal)])
    end
    
    % generate new Gamma for this annealing step
    if anneal_idx == 1
        Gamma = generate_feasible_Gamma(K,T,options.Gamma0);
    else
        Gamma = generate_feasible_Gamma(K,T,[]);
    end
    gamma_vec = reshape(Gamma',T*K,1); % vectorized form of gamma used in QP
    
    Lit = [];
    L = Inf;
    it = 0; % iteration counter
    while it < maxit % practical stopping criteria is present after computing new L (see "break")
        
        % compute Theta
        if isempty(options.Theta_given)
            for k=1:K
                sum_gammak = sum(Gamma(k,:));
                if sum_gammak ~= 0 % maybe gamma_k = 0 ? (i.e. this cluster is empty)
                    Theta(:,k) = sum(bsxfun(@times,X,Gamma(k,:)),2)/sum_gammak;
                else
                    Theta(:,k) = zeros(n,1);
                end
            end
        else
            Theta = options.Theta_given;
        end
        
        % compute new gamma
        % prepare new linear term in QP,
        % i.e. compute new residuum based on new theta
        g = 2*sum((kron(ones(1,K),X) - kron(Theta,ones(1,T))).^2,1)';
        
        % rescale linear term
        g = g/(T*n);
        
        % solve QP problem
        [gamma_vec,itQP] = spgqp(options.epssqr*HG, -g, gamma_vec, K, options.epssqr*normH, options.qp_eps, options.qp_maxit);
        
        Gamma = reshape(gamma_vec,T,K)';
        
        % compute new function value
        Lold = L; % store old function value
        Llin = dot(g,gamma_vec);
        Lquad = dot(HG*gamma_vec,gamma_vec);
        L = Llin + options.epssqr*Lquad;
        deltaL = Lold - L; % the progress of objective function, Lold > L (?)
        
        % display some informations, so we know that something is happening
        if options.debug
            disp([num2str(it) '. it: L = ' num2str(L) ', deltaL = ' num2str(deltaL) ', itQP = ' num2str(itQP)])
        end
        
        % stopping (breaking) criteria
        % based on sufficient decrease of objective funtion
        if abs(deltaL) < options.eps
            break; % stop outer "while" cycle
        end
        
        it = it + 1; % increase number of iterations
        
        if options.debug
            Lit(it) = L; % for postprocessing
        end
    end
    
    % check if the solution is better than previous one
    if options.nanneal > 1
        if L < L_best
            Theta_best = Theta;
            Gamma_best = Gamma;
            L_best = L;
            Llin_best = Llin;
            Lquad_best = Lquad;
            Lit_best = Lit;
        end
        
        if anneal_idx == options.nanneal
            % last annealing - return best values
            Theta = Theta_best;
            Gamma = Gamma_best;
            L = L_best;
            Llin = Llin_best;
            Lquad = Lquad_best;
            Lit = Lit_best;
        end
    end
end

end

