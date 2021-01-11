function E = smoE(alpha,X,y,b,K,index)
    E = 0;
    for i = 1:length(y)
        E = E + alpha(i)*y(i)*K(X(index,:),X(i,:));
    end
    E = E + b - y(index);
end