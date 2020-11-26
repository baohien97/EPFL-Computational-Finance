function xx = SimSDEJacobi(N_sim, N_time, T, X0, V0, kappa, sigma, theta,...
        rho, r, v_min, v_max)
    % We simulate the SDE for the Jacobi model by Euler discretization
    dt = T/N_time;
    xprev = X0*ones(N_sim,1);
    vprev = V0*ones(N_sim,1);
    t = 0;
    while t < T
        Z1 = randn(N_sim, 1);
        Z2 = randn(N_sim, 1);
        vnext = vprev + kappa*(theta - vprev)*dt + sigma*...
            sqrt(max(0,Q(vprev, v_min, v_max)))*sqrt(dt).*Z1;
        xnext = xprev + (r - vprev/2)*dt + rho*...
            sqrt(max(0,Q(vprev, v_min, v_max)))*sqrt(dt).*Z1 + ...
            sqrt(max(0,vprev - (rho^2)*Q(vprev, v_min, v_max)))*...
            sqrt(dt).*Z2;
        vprev = vnext;
        xprev = xnext;
        t = t + dt;
    end
    xx = xnext;
end

%%
function [q] = Q(v, v_min, v_max)
    diff_v_sq = (sqrt(v_max) - sqrt(v_min))^2;
    q = ((v - v_min).*(v_max - v))/diff_v_sq;
end
