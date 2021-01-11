function acc = binAcc(X,y,alpha,K,b)
    correct_count = 0;
    for i = 1:length(y)
        correct_count = correct_count + (y(i) == sign(value(X,y,alpha,K,b,X(i,:))));
    end
    acc = correct_count / length(y);
end