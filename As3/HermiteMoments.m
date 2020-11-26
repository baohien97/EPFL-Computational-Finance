%% main fn: Hermite moments
function [moments] = HermiteMoments(N, V0, X0, T, kappa, sigma, theta,...
        rho, r, mu_w, sigma_w, v_min, v_max)
    G = GenJacobi(N, kappa, sigma, theta, rho, r, sigma_w, v_min, v_max);
    
    % basic vectors B
    X0_tilde = (X0 - mu_w)/sigma_w;
    B = zeros(length(G), 1);
    for m=0:N
        for n=0:N-m
            ColInd = Ind(m,n);
            B(ColInd,1) = V0^m*HermitePol(n,X0_tilde);
        end
    end
    moments = B'*expm(G*T); 
end

%% Hermite polynomials (probabilistic)
function [H] = HermitePol(n, x)
    % probabilist hermite polynomial
    % x should be x_tilde = (x0-mu_w)/sigma_w
    H = (1/sqrt(factorial(n)))*2^(-0.5*n)*...
        hermiteH(n, x/sqrt(2));
end

%% G-matrix of Jacobi proc
function [G] = GenJacobi(N, kappa, sigma, theta, rho, r, sigma_w, v_min,...
        v_max)
    M = (N+2)*(N+1)/2;
    G = zeros(M, M);
    diff_v_sq = (sqrt(v_max) - sqrt(v_min))^2;
    for m = 0:N
        for n = 0:N-m
            ColInd = Ind(m,n);
            if m > 1
                G(Ind(m-2,n), ColInd) = -(sigma^2*m*(m-1)*v_max*v_min)/...
                    (2*diff_v_sq);
            end
            if m > 0 && n > 0
                G(Ind(m-1,n-1), ColInd) = -(sigma*rho*m*...
                    sqrt(n)*v_max*v_min)/(sigma_w*diff_v_sq);
            end
            if m > 0
                G(Ind(m-1,n), ColInd) = kappa*theta*m +...
                    sigma^2*m*(m-1)*(v_max+v_min)/(2*diff_v_sq);
            end
            if n > 0
                G(Ind(m,n-1), ColInd) = r*sqrt(n)/sigma_w +...
                    sigma*rho*m*sqrt(n)*(v_max + v_min)/(sigma_w*diff_v_sq);
            end
            if n > 1
                G(Ind(m+1,n-2), ColInd) = sqrt(n*(n-1))/(2*sigma_w^2);
            end
            if n > 0
                G(Ind(m+1,n-1), ColInd) = -sqrt(n)/(2*sigma_w) -...
                    sigma*rho*m*sqrt(n)/(sigma_w*diff_v_sq);
            end
            G(ColInd, ColInd) = -kappa*m - sigma^2*m*(m-1)/(2*diff_v_sq);
        end
    end
end

%% this is pi
function [ Index ] = Ind( m, n )
%Input:  m,n - powers of v and (x-muw)/sigmaw of basis elements.
% Output: Index - position of the corresponding basis element in Hn.
    Index = (m + n + 1)*(m + n) / 2 + n + 1;
end