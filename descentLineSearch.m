function [xMin, fMin, nIter, info] = descentLineSearch(F, descent, ls, alpha0, x0, tol, maxIter)
% DESCENTLINESEARCH Wrapper function executing  descent with line search
% [xMin, fMin, nIter, info] = descentLineSearch(F, descent, ls, alpha0, x0, tol, maxIter) 
%
% INPUTS
% F: structure with fields
%   - f: function handler
%   - df: gradient handler
%   - d2f: Hessian handler
% descent: specifies descent direction {'steepest', 'newton'}
% ls: function handle for computing the step length
% alpha0: initial step length 
% rho: in (0,1) backtraking step length reduction factor
% c1: constant in sufficient decrease condition f(x_k + alpha_k*p_k) > f_k + c1*alpha_k*(df_k')*p_k)
%     Typically chosen small, (default 1e-4). 
% x0: initial iterate
% tol: stopping condition on minimal allowed step
%      norm(x_k - x_k_1)/norm(x_k) < tol;
% maxIter: maximum number of iterations
%
% OUTPUTS
% xMin, fMin: minimum and value of f at the minimum
% nIter: number of iterations 
% info: structure with information about the iteration 
%   - xs: iterate history 
%   - alphas: step lengths history 
%
% Copyright (C) 2017  Marta M. Betcke, Kiko Rullan

% Parameters
% Stopping condition {'step', 'grad'}
stopType = 'step';

% Initialization
nIter = 0;
x_k = x0;
info.xs = x0;
info.alphas = alpha0;
stopCond = false; 

% Loop until convergence or maximum number of iterations
while (~stopCond && nIter <= maxIter)
    
    % ====================== YOUR CODE HERE ===================================
    % Instructions: x_k contains the current iteration point    - used in the stopping condition
    %               x_k_1 contains the previous iteration point - used in the stopping condition
    
    x_k_1 = x_k; % Set previous x_k
    if descent == "newton"
        p_k = -(F.d2f(x_k) \ F.df(x_k)); % Update x_k according to Newtons Method
        alpha0 = ls(x_k,p_k,alpha0); % Update step size
        x_k = x_k + alpha0 * p_k; % Update x_k in direction of gradient with step size
    else
        p_k = -F.df(x_k); % Get direction of gradient
        alpha0 = ls(x_k,p_k,alpha0); % Update step size
        x_k = x_k + alpha0 * p_k; % Update x_k in direction of gradient with step size
    end;
    
    info.xs = [info.xs;x_k]; % Append to list of xs
    info.alphas = [info.alphas;alpha0]; % Append to list of alphas
    nIter = nIter + 1; % Increment iterations

    % =========================================================================    
    switch stopType
      case 'step' 
        % Compute relative step length
        normStep = norm(x_k - x_k_1)/norm(x_k_1);
        stopCond = (normStep < tol);
      case 'grad'
        stopCond = (norm(F.df(x_k), 'inf') < tol*(1 + tol*abs(F.f(x_k))));
    end
    
end

% Assign output values 
xMin = x_k;
fMin = F.f(x_k); 

end