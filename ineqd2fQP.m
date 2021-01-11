function ineq = ineqd2fQP(alpha,C)
    ineq = zeros(length(alpha),length(alpha));
    for i = 1:length(alpha)
        for j = 1:length(alpha)
            if (i == j) && ((alpha(i) < 0) || (alpha(i) > C))
                ineq(i,j) = 2;
            end
        end
    end
end