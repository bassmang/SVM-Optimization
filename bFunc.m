function b = bFunc(alpha,X,y,C,K)
    b = 0;
    S = find((alpha > 0) & (alpha < C))';
    for s = S
        m_sum = 0;
        for m = S
            m_sum = m_sum + alpha(m)*y(m)*K(X(m,:),X(s,:));
        end
        b = b + y(s) - m_sum;
    end
    b = b / length(S);
end