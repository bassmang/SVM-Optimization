function [alpha,b,info] = svmQP(X,species,class,K,C,mu,muMult,t,tMult,tol,maxIter,getInfo)
    accs = [];
    y = getY(species,class);
    alpha0 = ones(length(y),1);
    H = zeros(length(y),length(y));
    for i = 1:length(y)
        for j = 1:length(y)
            H(i,j) = y(i)*y(j)*K(X(i,:),X(j,:));
        end
    end

    F.f = @(alpha) alpha'*H*alpha - sum(alpha);
    F.df = @(alpha) (H' + H)*alpha - 1;
    F.d2f = @(alpha) (H' + H);

    % Eq constraints
    eq.f = @(alpha) sum(y.*alpha')^2;
    eq.df = @(alpha) 2*sum(y.*alpha')*y';
    eq.d2f = @(alpha) 2*y'*y;

    % Combine constraints
    phi.f = @(alpha) eq.f(alpha) + ineqfQP(alpha,C);
    phi.df = @(alpha) eq.df(alpha) + ineqdfQP(alpha,C);
    phi.d2f = @(alpha) eq.d2f(alpha) + ineqd2fQP(alpha,C);

    % Set alpha and b
    [alpha, fMin, t, nIter, info] = quadraticPenalty(F, phi, alpha0, t, tMult, mu, muMult, tol, maxIter, X, y, C, K, getInfo);
    
    % Append Accuracy
    b = bFunc(alpha,X,y,C,K);
end