% part a: Monte-Carlo prices of up-and-out barrier over [0,T]

% r
% sigma
% N_time
% N_sm
% T 
% s
% K
% b

function [P_MC] = MCpriceBarrierUODM(r, sig, N_time, N_sim, T, s, K, b)
    %tk = T.*(0:N_time)/N_time;
    dt = T/N_time;
    S_half = s*ones(N_sim, 1);
    
    % half time
    for i=1:N_time/2 
        S_half = S_half.*(1 + r*dt + sig*sqrt(dt)*randn(N_sim,1));
    end
    
    % rest of the time
    S_end = S_half;
    for i=N_time/2+1:N_time
        S_end = S_end.*(1 + r*dt + sig*sqrt(dt)*randn(N_sim,1)); 
    end
    
    % prices at maturity
    P_T = max(0, S_end - K).*(S_half < b).*(S_end < b);
    
    % prices today
    P_MC = exp(-r*T)*mean(P_T);
    