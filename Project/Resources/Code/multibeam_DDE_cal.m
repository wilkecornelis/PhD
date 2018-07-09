%% multi-beam calibration of DDEs
%
% This script demonstrates the advantages of combining DD calibration
% solutions from multiple beams to improve the solutions for DDEs under the
% assumption that the DDE model can be split in a factor common to all
% beams and a global factor. This will be demonstrated using a simple model
% in which both factors can be described by low order polynomials.
%
% As starting point for the optimisation, we will assume that the apparent
% source fluxes are available / measured and that the noise is spatially
% white over the FoV.
%
% SJW, 29 January 2018

%% start with a clean workspace
clear
close all

%% measurement description
N_b = 9;        % number of beams (odd to ensure symmetric placement)
dl_b = 0.1;     % beam separation (assumed identical to HPBW)
dl = 1e-3;      % resolution for plots
Nrun = 1e3;     % number of runs in MC sim
Niter = 5;      % number of iterations
SNR = 10;       % instantaneous SNR of sources based on true power

%% define DD model per beam
p_b = polyfit([0, -dl_b/2, dl_b/2], [1, 0.5, 0.5], 2);
%l_b = -dl_b/2:dl:dl_b/2;
%g_b = polyval(l_b, p_b);

%% define global DD model
p_e = polyfit([0, -1, 1], [1, 0, 0], 2);
%l_e = -1:dl:1;
%g_e = polyval(p_e, l);

%% combine local and global DD model
l0 = -(N_b-1)/2 * dl_b:dl_b:(N_b-1)/2 * dl_b;    % beam centers
l_tot = ((l0(1) - dl_b/2):dl:(l0(end) + dl_b/2)).';
l_b = NaN(length(l_tot), 1);
%g_b = NaN(length(l_tot), 1);
for idx = 1:N_b
    l = (l0(idx) - dl_b/2):dl:(l0(idx) + dl_b/2);
    [~, startidx] = min(abs(l_tot - l(1)));
    [~, stopidx] = min(abs(l_tot -l(end)));
    l_b(startidx:stopidx) = l - l0(idx);
end
g_b = polyval(p_b, l_b);
g_e = polyval(p_e, l_tot);
g_tot = g_e .* g_b;
% figure
% plot(l_tot, g_b, 'k:', ...
%      l_tot, g_e, 'k-.', ...
%      l_tot, g_tot, 'k-');
% set(gca, 'FontSize', 16);
% xlabel('l (direction cosine)');
% ylabel('gain (normalised to unity)');

%% perform MC sim
lidx = 1:10:length(l_tot);
%lidx = randperm(length(l_tot), 3 * N_b);        % grid indices of sources
sigma0 = SNR;               % true source power
g_global = NaN(length(l_tot), Nrun);
g_local2 = NaN(length(l_tot), Nrun);
g_local4 = NaN(length(l_tot), Nrun);
for run = 1:Nrun
    sigma = g_tot(lidx) .* (sigma0 + randn(length(lidx), 1));
    % try iterative multi-beam solution
    g_e_hat = ones(length(l_tot), 1);
    g_b_hat = ones(length(l_tot), 1);
    for iter = 1:Niter
        % remove (estimated) global DD model
        sigma_hat = sigma ./ g_e_hat(lidx);
        % update local DD model
        p_b_hat = polyfit(l_b(lidx), sigma_hat, 2);
        g_b_hat = polyval(p_b_hat, l_b);
        % remove (estimated) local DD model
        sigma_hat = sigma ./ g_b_hat(lidx);
        % update global DD model
        p_e_hat = polyfit(l_tot(lidx), sigma_hat, 2);
        g_e_hat = polyval(p_e_hat, l_tot);
    end
    g_global(:, run) = g_e_hat .* g_b_hat;
    % try local solution (per beam solution)
    for idx = 1:N_b
        l = (l0(idx) - dl_b/2):dl:(l0(idx) + dl_b/2);
        [~, startidx] = min(abs(l_tot - l(1)));
        [~, stopidx] = min(abs(l_tot - l(end)));
        srcidx = find((l_tot(lidx) >= l(1)) & (l_tot(lidx) <= l(end)));
        p_local2 = polyfit(l_tot(lidx(srcidx)) - l0(idx), sigma(srcidx), 2);
        g_local2(startidx:stopidx, run) = polyval(p_local2, l_tot(startidx:stopidx) - l0(idx)); 
        p_local4 = polyfit(l_tot(lidx(srcidx)) - l0(idx), sigma(srcidx), 4);
        g_local4(startidx:stopidx, run) = polyval(p_local4, l_tot(startidx:stopidx) - l0(idx)); 
    end
end

%% plot statistics
g_global = g_global / SNR;
g_local2 = g_local2 / SNR;
g_local4 = g_local4 / SNR;
std_global = std(g_global.');
bias_global = mean((g_global - g_tot * ones(1, Nrun)).');
std_local2 = std(g_local2.');
std_local4 = std(g_local4.');
bias_local2 = mean((g_local2 - g_tot * ones(1, Nrun)).');
bias_local4 = mean((g_local4 - g_tot * ones(1, Nrun)).');
figure
plot(l_tot, std_local2, 'k:', ...
     l_tot, std_local4, 'k-.', ...
     l_tot, std_global, 'k-');
set(gca, 'FontSize', 16);
xlabel('l (direction cosine)');
ylabel('std (gains normalised to unity)');
figure
plot(l_tot, bias_local2, 'k:', ...
     l_tot, bias_local4, 'k-.', ...
     l_tot, bias_global, 'k-');
set(gca, 'FontSize', 16);
xlabel('l (direction cosine)');
ylabel('bias (gains normalised to unity)');
