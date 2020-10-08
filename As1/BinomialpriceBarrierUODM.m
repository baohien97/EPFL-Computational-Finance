% part a: multi-period binomial model over [0,T] with the given params

% u: up factor
% d: down factor
% N: number of periods
% T: maturity
% s: S0 - initial price of underlying
% K: strike of underlying
% barrier

function  [P_bin] = BinomialpriceBarrierUODM(r, d, u, N, T, s, K, b)
    
    % underlying tree
    S_tree = nan(N+1,N+1);
    S_tree(1,1) = s;
    for i=2:N+1
        S_tree(1:i-1,i) = S_tree(1:i-1,i-1)*u; % up nodes
        S_tree(i,i) = S_tree(i-1,i-1)*d; % down nodes
    end
    
    % option tree
    r_tilde = r*T/N;
    q_u = ((1+r_tilde)-d)/(u-d);
    q_d = (u - (1+r_tilde))/(u-d);
    S_tree(isnan(S_tree)) = 0;
    P_tree = zeros(size(S_tree));
    
    % monitoring date N
    for i=1:N+1
        if (S_tree(i, end) < b)
            P_tree(i, end) = max(0, S_tree(i, end) - K);
        end
    end
    
    for i=N:-1:N/2+1
        P_tree(1:i,i) = (1/(1+r_tilde))*(q_u*P_tree(1:i,i+1) + ...
            q_d*P_tree(2:i+1,i+1));
    end
    
    % monitoring date N/2
    for i=1:N/2+1
        if (S_tree(i, N/2 + 1) > b)
            P_tree(i, N/2 +1) = 0;
        end
    end
    
    % pricing option from 0 to N/2
    for i=N/2:-1:1
        %disp(N)
        P_tree(1:i,i) = (1/(1+r_tilde))*(q_u*P_tree(1:i,i+1) + ...
            q_d*P_tree(2:i+1,i+1));
    end
    
    %P_tree(isnan(P_tree)) = 0;
    P_bin = P_tree(1,1);
