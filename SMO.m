function [alpha,b,info] = SMO(X,species,class,K,C,tol,max_passes,maxIter,getInfo)

    y = getY(species,class);
    iters = 0;
    m = length(y);
    alpha = zeros(m,1);
    b = 0;
    passes = 0;
    info.accs = [];
    info.xnorms = [];
    info.fnorms = [];

    while (passes < max_passes && iters < maxIter)
        num_changed_alphas = 0;
        alpha_old = alpha;
        for i = 1:m
            E_i = smoE(alpha,X,y,b,K,i);
            if (((y(i)*E_i < -tol) && (alpha(i) < C)) || ((y(i)*E_i > tol) && (alpha(i) > 0)))
                all_inds = 1:m; % Get all indices
                all_inds(i) = []; % Remove current index
                j = randsample(all_inds,1); % Get random j, j!= i
                alpha_i_old = alpha(i);
                alpha_j_old = alpha(j);
                if y(i) == y(j)
                    L = max(0,alpha(i) + alpha(j) - C);
                    H = min(C,alpha(i) + alpha(j));
                else
                    L = max(0, alpha(j) - alpha(i));
                    H = min(C, C + alpha(j) - alpha(i));
                end
                if (H == L)
                    continue;
                end
                eta = 2*K(X(i,:),X(j,:)) - K(X(i,:),X(i,:)) - K(X(j,:),X(j,:));
                if (eta >= 0)
                    continue;
                end
                E_j = smoE(alpha,X,y,b,K,j);
                alpha_j_new = alpha(j) - (y(j)*(E_i - E_j) / eta);
                if (alpha_j_new > H)
                    alpha(j) = H;
                elseif (alpha_j_new < L)
                    alpha(j) = L;
                else
                    alpha(j) = alpha_j_new;
                end
                if (abs(alpha(j) - alpha_j_old) < 1E-5)
                    continue;
                end
                alpha(i) = alpha(i) + y(i)*y(j)*(alpha_j_old - alpha(j));
                b1 = b - E_i - y(i)*(alpha(i) - alpha_i_old)*K(X(i,:),X(i,:)) - y(j)*(alpha(j) - alpha_j_old)*K(X(i,:),X(j,:));
                b2 = b - E_j - y(i)*(alpha(i) - alpha_i_old)*K(X(i,:),X(j,:)) - y(j)*(alpha(j) - alpha_j_old)*K(X(j,:),X(j,:));
                if (alpha(i) > 0) && (alpha(i) < C)
                    b = b1;
                elseif (alpha(j) > 0) && (alpha(j) < C)
                    b = b2;
                else
                    b = (b1 + b2) / 2;
                end
                num_changed_alphas = num_changed_alphas + 1;
            end
        end
        if (num_changed_alphas == 0)
            passes = passes + 1;
        else
            passes = 0;
        end
        
        % Get info if flag set
        if getInfo
            info.xnorms = [info.xnorms norm(alpha - alpha_old)];
            info.fnorms = [info.fnorms norm(SMOObj(X,y,alpha,K) - SMOObj(X,y,alpha_old,K))];
            info.accs = [info.accs binAcc(X,y,alpha,K,b)];
        end
        
        % Update iteration
        iters = iters + 1;
        
    end
end