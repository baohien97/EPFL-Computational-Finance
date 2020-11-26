X0 = 0;
V0 = 0.04;
kappa = 0.5;
theta = 0.04;
sigma = 1;
r = 0;
rho = -0.5;
T = 1/12;
v_min = 10^(-4);
v_max = 0.08;
N = 10;

% for heston
L = 100;
alpha = 1;

ks = linspace(-0.1,0.1,50);
imp_vols_jacobi = zeros(length(ks),1);
imp_vols_heston = zeros(length(ks),1);

for i=1:length(ks)
    jacobi_price = PriceApprox(N, V0, X0, T, ks(i), kappa, sigma, theta,...
        rho, r, v_min, v_max); 
    exp_k = exp(ks(i));
    %disp(ks(i));
    %disp(jacobi_price);
    imp_vols_jacobi(i) = blsimpv(exp(X0), exp_k, r, T, jacobi_price);
    heston_price = CallHestonFourier(exp(X0), exp_k, T, r, kappa, theta,...
        sigma, rho, V0, alpha, L); 
    imp_vols_heston(i) = blsimpv(exp(X0), exp_k, r, T, heston_price);
end


%% plotting 
figure
plot(ks, imp_vols_jacobi, "b", ks, imp_vols_heston, "r");
title('Volatility smile')
xlabel('k=log(K)')
ylabel('implied vol')
legend("Jacobi", "Heston")

%% Heston Fourier (ex 2 Sheet 6)
function [P] = CallHestonFourier(S,K,T,r,kappa,theta,sigma,rho,V,alpha,L)
    b=@(nu)(kappa-1i*rho*sigma.*nu);
    gamma=@(nu)(sqrt(sigma^2*(nu.^2+1i.*nu)+b(nu).^2));
    a=@(nu)(b(nu)./gamma(nu)).*sinh(T*0.5.*gamma(nu));
    c=@(nu)(gamma(nu).*coth(0.5*T.*gamma(nu))+b(nu));
    d=@(nu)(kappa*theta*T.*b(nu)/sigma^2);
    f=@(nu)(1i*(log(S)+r*T).*nu+d(nu));
    g=@(nu)(cosh(T*0.5.*gamma(nu))+a(nu)).^(2*kappa*theta/sigma^2);
    h=@(nu)(-(nu.^2+1i.*nu)*V./c(nu));
    phi=@(nu)(exp(f(nu)).*exp(h(nu))./g(nu));
    
    integrand=@(nu)(real((phi(nu-1i*(alpha+1))./...
        (alpha^2+alpha-nu.^2+1i*(2*alpha+1)*nu)).*exp(-1i*log(K).*nu))); ...
        % original integrand

    P=(exp(-r*T-alpha*log(K))/pi)*integral(integrand,0, L);% Price
end