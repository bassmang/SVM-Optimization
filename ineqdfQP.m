function ineq = ineqdfQP(alpha,C)
    ineq = zeros(length(alpha),1);
    for i = 1:length(alpha)
        if alpha(i) < 0
            ineq(i) = 2*alpha(i);
        elseif alpha(i) > C
            ineq(i) = 2*(C - alpha(i))*(-1);
        end
    end
end