%% Given params
S0 = 1;
lambda = 0.4;
sigma = 0.15;
alpha = -.5;
beta = .4;
T = .5;
a_min  = 0.7;
a_max = 1.3;
L = 50;
eta = -1;

a = linspace(a_min, a_max, 100);

fourrier_prices = zeros(length(a),1);

x0 = log(S0);
start_f = tic;
for i=1:length(a)
    fourrier_prices(i) = EuropeanDigital(x0, T, eta, a(i), L,...
        sigma, alpha, beta, lambda);
end
end_f = toc;

orders = 2:1:30;
mae_n = zeros(length(orders),1);
cheb_prices = zeros(length(a),length(orders));
cheb_ex_times = length(orders);

for i=1:length(orders)
    tic
    f = @(a) EuropeanDigital(x0, T, eta, a, L, sigma, alpha,...
        beta, lambda);
    cheb_prices(:,i) = ChebInterpol(f, a, orders(i), a_min, a_max);
    cheb_ex_times(:,i) = toc;
    mae_n(i) = max(abs(fourrier_prices - cheb_prices(:, i)));
end

%% plots
figure(1)
plot(a, fourrier_prices, "r-", a, cheb_prices(:,end), "b*");
xlabel("a");
ylabel("price");
title("Prices")
legend('Fourrier Prices', 'Chebyshev-interpolated prices');
saveas(gcf, "Prices.png");

figure(2)
plot(orders, mae_n);
xlabel("order n");
ylabel("maximal absolute error")
title("Maximal Absolute Errors")
saveas(gcf, "mae.png");

figure(3)
plot(orders, cheb_ex_times, "r-");
xlabel("order n");
ylabel("seconds");
title("Execution time (Chebyshev)")
saveas(gcf, "exec_time.png");

