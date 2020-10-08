%% Params
s = 1;
r = .1;
T = .5;
K = .9;
sig = .1;
b = 1.3;

%% MC prices
N_time = 100;
N_sim = 1.0e6;
MCprice = MCpriceBarrierUODM(r, sig, N_time, N_sim, T, s, K ,b);

%% Binomial prices
N = (2:2:200);
binprices = zeros(length(N),1);
for i=1:length(N)   
    u = 1 + r*T/N(i) + sig*sqrt(T)/sqrt(N(i));
    d = 1 + r*T/N(i) - sig*sqrt(T)/sqrt(N(i));
    binprices(i) = BinomialpriceBarrierUODM(r, d, u, N(i), T, s, K, b);
end

%% Plots
figure
plot(N, MCprice*ones(1,length(N)), "r")
hold on
plot(N, binprices, "b");
hold off

legend({"Monte-Carlo", "Binomial"}, "Location", "bestoutside")
xlabel("N")
ylabel("Option Price")
title("Prices of Up-and-Out option")