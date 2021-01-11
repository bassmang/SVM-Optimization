function obj = SMOObj(X,y,alpha,K)
	lag_sum = 0;
    for i = 1:length(y)
        for j = 1:length(y)
            lag_sum = lag_sum + y(i)*y(j)*alpha(i)*alpha(j)*K(X(i,:),X(j,:));
        end
    end
    obj = .5*lag_sum - sum(alpha);
end