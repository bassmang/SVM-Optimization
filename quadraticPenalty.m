function [xMin, fMin, t, nIter, info] = quadraticPenalty(F, phi, x0, t, tMult, mu, muMult, tol, maxIter, X, y, C, K, getInfo)
nIter = 0;
stopCond = false;
x_k = x0;
info.xnorms = [];
info.fnorms = [];
info.accs = [];

% Parameters for centering step
alpha0 = 1; 
opts.c1 = 1e-4;
opts.c2 = 0.9;
opts.rho = 0.5;
maxIterNewton = 100;
% Loop 

while (~stopCond && nIter < maxIter)
    % Create function handler for centering step (needs to be redifined at each step because of changing "mu")
    G.f = @(x) F.f(x) + (mu/2)*phi.f(x);
    G.df = @(x) F.df(x) + (mu/2)*phi.df(x);
    G.d2f = @(x) F.d2f(x) + (mu/2)*phi.d2f(x);
    % Line search function (needs to be redefined at each step because of changing G) 
    lsFun = @(x_k, p_k, alpha0) backtracking(G, x_k, p_k, alpha0, opts);
    % Centering step
    x_k_1 = x_k;
    [x_k, f_k, nIterLS, infoIter] = descentLineSearch(G, 'gradient', lsFun, alpha0, x_k, t, maxIterNewton);   
    % Check stopping condition (m/t). Assumes the tolerance has been scaled with 1/m
    if norm(x_k - x_k_1) < tol; stopCond = true; end

    % Increase mu at each iteration
    mu = mu*muMult;
    
    % Decrease t at each iteration
    t = t*tMult;

    % Include info if flag is set
    if getInfo
        b = bFunc(x_k,X,y,C,K);
        info.accs = [info.accs binAcc(X,y,x_k,K,b)];
        info.xnorms = [info.xnorms norm(x_k - x_k_1)];
        info.fnorms = [info.fnorms norm(G.f(x_k) - G.f(x_k_1))];
    end
    
    % Increment number of iterations
    nIter = nIter + 1;
    
end

% Assign values
xMin = x_k;
fMin = F.f(x_k);
end

