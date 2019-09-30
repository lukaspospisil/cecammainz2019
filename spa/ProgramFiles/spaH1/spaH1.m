function [ S, Gamma, L, Llin, Lquad, Lit ] = spaH1( X,Hreg,options )

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

% if box constraints are given only for one dimension, then copy it to all
if ~isempty(isempty(options.Sbox))
    if and(n > 1, size(options.Sbox,1) == 1)
        options.Sbox = kron(options.Sbox, ones(n,1));
    end
end

HG = kron(speye(K),Hreg);

% here will be stored solution of model parameters - one for each cluster
S = zeros(n,K);

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
    
    L = Inf;
    Lit = []; % function value through iterations
    it = 0; % iteration counter
    while it < maxit % practical stopping criteria is present after computing new L (see "break")
        
        % compute Theta
        if isempty(options.Theta_given)
            A = Gamma*Gamma' + (2*options.Sreg/(n*K*(K-1)))*(K*eye(K) - ones(K,K));
            B = Gamma*X';
            
            if isempty(options.Sbox)
                % solve system(s) of linear equations
                S = solve_S(A,B,options.Sreg > 1e-8);
            else
                % solve QP with box
                S = solve_Sbox(A,B,options.Sbox);
            end
        else
            S = options.Theta_given;
        end
        
        % compute new gamma
        HS=S'*S; HS=(HS+HS');
        H = 2*options.epssqr*HG + kron(HS,speye(T))/(T*n);
        CTX = 2*reshape((X'*S),K*T,1)/(T*n);
        
        normH = gershgorin(H); %eigs(H,1);
        
        % solve qp
        [gamma_vec,itQP] = spgqp(H, CTX, gamma_vec, K, normH, options.qp_eps, options.qp_maxit);
        
        Gamma = reshape(gamma_vec,T,K)';
        
        % compute new function value
        Lold = L; % store old function value
        
        Llin = norm(X - S*Gamma,'fro')^2/(T*n);
        Lquad = gamma_vec'*HG*gamma_vec;
        L = Llin + options.epssqr*Lquad;
        deltaL = Lold - L; % the progress of objective function, Lold > L (?)
        
        it = it + 1; % increase number of iterations of subspace algorithm
        
        % display some informations, so we know that something is happening
        if options.debug
            disp([num2str(it) '. it: L = ' num2str(L) ', deltaL = ' num2str(deltaL) ', itQP = ' num2str(itQP)])
            Lit(it) = L;
        end
        
        % stopping (breaking) criteria
        % based on sufficient decrease of objective funtion
        if abs(deltaL) < options.eps
            break; % stop outer "while" cycle
        end
        
    end
    
    % check if the solution is better than previous one
    if options.nanneal > 1
        if L < L_best
            S_best = S;
            Gamma_best = Gamma;
            L_best = L;
            Llin_best = Llin;
            Lquad_best = Lquad;
            Lit_best = Lit;
        end
        
        if anneal_idx == options.nanneal
            % last annealing - return best values
            S = S_best;
            Gamma = Gamma_best;
            L = L_best;
            Llin = Llin_best;
            Lquad = Lquad_best;
            Lit = Lit_best;
        end
    end
    
end

end

function Theta = solve_S(A,B,flag_isreg)
% system is regular with regularization
if flag_isreg
    Theta = ( A\B )';
else
    Theta = ( pinv(A)*B )';
end
end

function Theta = solve_Sbox(A,B,box_constraints)

[K,n] = size(B);

options = optimoptions('quadprog','Display','none');

Theta = zeros(n,K);
for i = 1:n
    lb = box_constraints(i,1)*ones(K,1);
    ub = box_constraints(i,2)*ones(K,1);
    
    Theta(i,:) = quadprog(A,-B(:,i)',[],[],[],[],lb,ub,[],options);
end

end