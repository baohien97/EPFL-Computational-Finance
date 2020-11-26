X0 = 0;
V0 = 0.04;
kappa = 0.5;
theta = 0.04;
sigma = 1;
r = 0;
rho = -0.5;
T = 1/12;
v_min = 1e-5;
v_max = 0.08;
k = -0.1; % log strike


%% polynomial expansion
N = 50;
price_poly = PriceApprox(N,V0,X0,T,k,kappa,sigma,theta,rho,r,...
    v_min,v_max);


%% Monte-carlo
exp_k = exp(k);
N_sim = 10^6;
N_time = 100;

X_sim = SimSDEJacobi(N_sim,N_time,T,X0,V0,kappa,sigma,theta,rho,r,...
    v_min,v_max); % log underlying at maturity
S_T_sim = exp(X_sim);
payoffs = max(0, S_T_sim - exp_k);

price_sim = mean(exp(-r*T)*payoffs);
