global T K S r C

%% Initial Params
theta = 0.04;
kappa = 1.50;
sigma = 0.30;
rho = -0.6;
V0 = 0.0441;

x0 = [theta, kappa, sigma, rho, V0]; % starting values
lb = [-Inf, 0, 0, -1, 0]; % lower bounds
ub = [Inf, Inf, Inf, 1, Inf]; % upper bounds

%% Other already-given variables
r = 0.015;
S = 1202.10;

data = load('Call_20050103.mat');
data = data.Call_20050103;
K = data(:, 2);
T = data(:, 3)/252;
C = data(:, 1);
% imp_vol = data(:, 4); => not needed in this assignment

%% Optimisation 
[x,residual,exitflag,output] = fminsearchcon(@(x) costfn(x), x0, lb, ub);
disp('Opt 2: [nu, kappa, sigma, rho, V0] = ')
disp(x)
disp('Num of iters: ')
disp(output.iterations)

%% Cost fn

function [cost] = costfn(x)
    global T K S r C
    cost = 0;
    for i=1:length(C)
        cost = cost + (Call_Heston(K(i), T(i), r, x(1), x(2), x(3), x(4),...
                S, x(5)) - C(i))^2;
    end
    cost = sqrt(mean(cost));
    %disp(cost)
end

