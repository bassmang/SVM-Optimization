function ineq = ineqfQP(alpha,C)
    ineq = 0;
    for i = 1:length(alpha)
        if alpha(i) < 0
            ineq = ineq + alpha(i)^2;
        elseif alpha(i) > C
            ineq = ineq + (C - alpha(i))^2;
        end
    end
end