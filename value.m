function v = value(X,y,alpha,K,b,z)
    v = 0;
    for i = 1:length(y)
        v = v + alpha(i)*y(i)*K(X(i,:),z);
    end
    v = v + b;
end