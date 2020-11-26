%% main fn
function [price] = PriceApprox(N, V0, X0, T, k, kappa, sigma, theta, rho,...
        r, v_min, v_max)

    % In order to determine muw and sigma_w, we first assume mu_w = 0 and
    % sigma_w = 1 and we apply the moment formula to compute E[X_T] and 
    % Var[X_T]
    sigma_w = 1;
    G = GenJacobi(2, kappa, sigma, theta, rho, r, sigma_w, v_min, v_max);
    H = [1, V0, X0, V0^2, V0*X0, X0^2];
    moments = H*expm(G*T);
    mu_w = moments(3);
    sigma_w = sqrt(moments(6) - mu_w^2);
    
    % Hermite moments
    moments = HermiteMoments(N,V0,X0,T,kappa,sigma,theta,rho,r,mu_w,...
        sigma_w,v_min,v_max);
    coefs = FourierCoefficients(N,k,mu_w,sigma_w,r,T);
    %disp(coefs);
    price_call = 0;
    for n=0:N
        price_call = price_call + moments(Ind(0,n))*coefs(n+1);
    end
    price = price_call;
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