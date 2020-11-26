%clear all;
clc
clear
%% params
sigma = 0.2;
r = 0.015;
K = 1;
T = 1;
S_min = 0;

trunc_err = 10^(-6);

% follow homework 7
S_max = K*exp(sigma^2*T/2-sigma*sqrt(T)*norminv(trunc_err/K)); 

theta = 0.5;
N_x = 60;
N_t = 100;

% forcing
rhs = @(S,t) zeros(length(S), length(t));

% bc right
bc_right = @(t) zeros(length(t),1);

% bc left
%bc_left = @(S,t) 0;

% initial cond
initial_cond = @(S) max(K - S, 0);

%% part b: compute V (original PDE)
[V,FD_grid,time_steps] = bs_timestepping(sigma,r,rhs,bc_right,...
    initial_cond,S_max,N_x,T,N_t,theta);
U_rhs = vega_rhs(V, FD_grid, sigma);

vega_bc_right = @(t) zeros(length(t),1);
vega_bc_left = @(S) zeros(length(S),1);

[Vega,FD_grid,~] = bs_vega_timestepping(sigma,r,U_rhs,vega_bc_right,...
    vega_bc_left,S_max,N_x,T,N_t,theta);

Vega_exact = blsvega(FD_grid(2:end),K,r,T,sigma);

% plot 
figure
plot(FD_grid(2:end),Vega(2:end,end),'r--',"LineWidth",1)
hold on
plot(FD_grid(2:end),Vega_exact,'b')
xlabel('S')
ylabel('Vega')
legend('Numerical','Exact')
saveas(gcf,"Ex2-1.png")

%% part c

N_xs = 10*2.^(0:5);
N_ts = 0.6.*N_xs;
errors = zeros(length(N_xs),1);
for i=1:length(N_xs)
        % time discretization parameters
    N_t = N_ts(i);
    N_x = N_xs(i);
    % solve for V
    [V,FD_grid] = bs_timestepping(sigma,r,rhs,bc_right,...
        initial_cond,S_max,N_x,T,N_t,theta);

    U_rhs = vega_rhs(V, FD_grid, sigma);

    % solve for vega
    [Vega,FD_grid,time_steps] = bs_vega_timestepping(sigma,r,U_rhs,...
        vega_bc_right,vega_bc_left,S_max,N_x,T,N_t,theta);

    % exact solution
    Vega_exact = blsvega(FD_grid(2:end)',K,r,T,sigma);
    
    % error
    errors(i) = max(abs(Vega(2:end,end)-Vega_exact));
    
end

P = polyfit(log(h),log(errors),1);
disp('order of convergence:' + string(P(1)))

% plot
h = S_max./N_xs;
figure(2)
plot(h, errors);
xlabel("h");
ylabel("error")
saveas(gcf,"Ex2-2.png")

figure(3)
loglog(h, errors);
xlabel("h");
ylabel("error")
saveas(gcf,"Ex2-2-loglog.png")

%% part d: central finite difference scheme
delta_sigma = [0.0005 0.001 0.01 0.02 0.03 0.05];
errors_d = zeros(length(N_xs), length(delta_sigma));
for j=1:length(delta_sigma)
    for i=1:length(N_xs)
        [V_up,FD_grid,~] = bs_timestepping(sigma+delta_sigma(j),r,rhs,...
            bc_right, initial_cond,S_max,N_xs(i),T,N_ts(i),theta);
        [V_down,FD_grid,~] = bs_timestepping(sigma-delta_sigma(j),r,...
            rhs,bc_right,initial_cond,S_max,N_xs(i),T,N_ts(i),theta);
        
        Vega = 0.5*(V_up(:,end) - V_down(:,end))/delta_sigma(j);
        errors_d(i,j) = max(abs(Vega(2:end) - blsvega(FD_grid(2:end)',...
            K,r,T,sigma)));
    end
end

% plot
figure(4)
plot(h, errors_d);
xlabel("h");
ylabel("error")
legend("ds = " + string(delta_sigma), 'Location', 'southeast')
saveas(gcf,"Ex2-3.png")

figure(5)
loglog(h, errors_d);
xlabel("h");
ylabel("error")
legend("ds = " + string(delta_sigma), 'Location', 'southeast')
saveas(gcf,"Ex2-3-loglog.png")
%% aux fn
function [rhs] = vega_rhs(U, FD_grid, sigma)
    [Nx, Nt] = size(U);
    h = FD_grid(end)/(Nx-1);
    
    rhs = zeros(Nx, Nt);
    for s=2:Nx-1
        % slide 27,28 lecture 7
        rhs(s,:) = sigma*FD_grid(s)^2/h^2 *(U(s+1,:) - 2*U(s,:) +...
            U(s-1,:));
    end
end