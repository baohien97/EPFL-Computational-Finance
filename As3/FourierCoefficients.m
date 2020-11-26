function [coefficients] = FourierCoefficients(N, k, mu_w, sigma_w, r, T)
    mu_tilde = (k-mu_w)/sigma_w;
    
    I0 = exp(sigma_w^2/2)*normcdf(sigma_w - mu_tilde);
    I = zeros(N+1, 1);
    I(1) = I0;
    
    % vector of coefficients
    f0 = exp(-r*T+mu_w)*I(1) - exp(-r*T+k)*normcdf(-mu_tilde);
    f = zeros(N+1, 1);
    f(1) = f0;
    
    for n=2:N+1
        % hermite poly takes order n-2 because index in matlab starts at 1
        I(n) = StandardHermitePol(n-2, mu_tilde)*exp(sigma_w*mu_tilde)*...
            normpdf(mu_tilde) + sigma_w*I(n-1);
        f(n) = exp(-r*T+mu_w)*(1/sqrt(factorial(n-1)))*sigma_w*I(n-1);
    end
    coefficients = f;
end

%% Standard Hermite polynomials (probabilist)
function [H] = StandardHermitePol(n, x)
    % probabilist hermite polynomial
    % x should be x_tilde = (x0-mu_w)/sigma_w
    H = 2^(-0.5*n)*hermiteH(n, x/sqrt(2));
end