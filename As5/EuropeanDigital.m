function [P] = EuropeanDigital(x0, T, eta, a, L, sigma, alpha,...
        beta, lambda)
     gamma = -0.5*sigma^2 - lambda*(exp(alpha + 0.5*beta^2) - 1);
 
     g_hat = @(z) 1i *exp(1i * log(a)*z) ./ z;
     P_hat = @(z) (exp(T*(1i*z*gamma...
         - 0.5*(sigma^2)*z.^2 + lambda*(exp(1i*z*alpha - 0.5*(beta^2)*z.^2)...
         - 1))));
 
     fun = @(xi) real(exp(1i*xi*x0).*g_hat(-xi-1i*eta).*P_hat(xi+1i*eta));
     
     % integration routine here
     P = exp(-eta*x0)/(2*pi)*integral(fun,0,L);
end


