%% Clear all things
clc; clear; close all; path(pathdef);
addpath('baselines/PrecondSPA/')
addpath('~/code/matlab/common/prob_tools')
addpath('utils')
addpath('~/code/matlab/common/PGD/')
addpath('~/code/matlab/common/prox_ops')
addpath('fw_core')
addpath('baselines/FGNSR/matlab')

M=50; N=100; K=20;


RAE_fw = [];
RAE_spa = [];
RAE_fg = [];

succ_fw = [];
succ_spa = [];
succ_fg = [];

succ2_fw = [];
succ2_spa = [];
succ2_fg = [];

MRSA_fw = [];
MRSA_spa = [];
MRSA_fg = [];

% list_eps = [0.01 0.02 0.05 0.08 0.1 0.2 0.3 0.4 0.5 0.8 1];
list_SNR = [5:1:30];

t0 = tic;
for exp_ind=1:numel(list_SNR)
    fprintf('%d/%d - Elapsed time: %f \n', exp_ind, numel(list_SNR), toc(t0));

    for trial_ind=1:5
        W = rand(M, K);
        W = W./sum(W, 1);

        H = zeros(K, N);
        H(:, 1:K) = eye(K);
        for col=K+1:N
            inds = randperm(K);
            tmp = dirichlet_rnd([50 50], 1);

            H(inds(1), col) = tmp(1);
            H(inds(2), col) = tmp(2);
        end
        % w_bar = mean(W, 2) + (1e-1)*randn(M, N-K);
        w_bar = mean(W, 2);

        Y = W*H;
        pure_pixel_set = 1:K;

        original_noise = zeros(M, N);
        original_noise(:, K+1:end) = Y(:, K+1:end) - w_bar;
        original_noise = original_noise + (1e-1)*randn(M, N);

        noise_energy = sqrt(norm(Y, 'fro')^2/10^(list_SNR(exp_ind)/10));
        noise = original_noise/norm(original_noise, 'fro')*noise_energy;
        % noise = original_noise/norm(original_noise, 'fro')*list_eps(j);
        X = Y + noise;
        SNR(exp_ind, trial_ind) = 10*log10(norm(Y, 'fro')^2/norm(noise, 'fro')^2);

        new_order = randperm(N);
        r_pure_pixel_set = [];
        for i=1:numel(pure_pixel_set)
            r_pure_pixel_set(end+1) = find(new_order == pure_pixel_set(i));
        end
        pure_pixel_set = r_pure_pixel_set;
        X = X(:, new_order);

        % ------------- SPA
        Lambda_hat = SPAselect(X, K);
        succ_spa(exp_ind, trial_ind) = all(sort(Lambda_hat(:)) == sort(pure_pixel_set(:)));
        succ2_spa(exp_ind, trial_ind) = numel(intersect(Lambda_hat, pure_pixel_set))/K;
        W_hat = X(:, Lambda_hat);
        H_hat = computeApproxError(X, Lambda_hat);
        RAE_spa(exp_ind, trial_ind) = 1-norm(X - W_hat*H_hat, 'fro')/norm(X, 'fro');
        MRSA_spa(exp_ind, trial_ind) = compute_MRSA(W_hat, W);


        % ------------- FRANK-WOLFE
        fw_tic = tic;
        options = struct;
        options.maxIters = 100;
        options.verbose = 0;
        options.lambda = -3;
        options.debug = 0;
        options.N = K;
        options.backend = 'mex';
        options.dualGapThreshold=0;
        options.epsilon = 1e-5;
        [C_hat, fw_tracking] = fw(X, options);
        [~, Lambda_hat] = maxk(vecnorm(C_hat, Inf, 2), K, 1);

        succ_fw(exp_ind, trial_ind) = all(sort(Lambda_hat(:)) == sort(pure_pixel_set(:)));
        succ2_fw(exp_ind, trial_ind) = numel(intersect(Lambda_hat, pure_pixel_set))/K;

        W_hat = X(:, Lambda_hat);
        H_hat = computeApproxError(X, Lambda_hat);
        RAE_fw(exp_ind, trial_ind) = 1-norm(X - W_hat*H_hat, 'fro')/norm(X, 'fro');
        MRSA_fw(exp_ind, trial_ind) = compute_MRSA(W_hat, W);



        [C_hat_fg, K_fg, ~, ~, ~, ~, t] = fgnsr(X, K, 'maxiter', 150, ...
            'debug', 0, 'verbose', false);
        succ_fg(exp_ind, trial_ind) = all(sort(K_fg(:)) == sort(pure_pixel_set(:)));
        succ2_fg(exp_ind, trial_ind) = numel(intersect(K_fg, pure_pixel_set))/K;

        W_hat = X(:, K_fg);
        H_hat = computeApproxError(X, K_fg);
        RAE_fg(exp_ind, trial_ind) = 1-norm(X - W_hat*H_hat, 'fro')/norm(X, 'fro');
        MRSA_fg(exp_ind, trial_ind) = compute_MRSA(W_hat, W);
    end
end

% save('./results/gillis_setting_SNR.mat')
% load('./results/gillis_setting_SNR.mat')
% figure();
% semilogx(list_eps, mean(RAE_spa, 2), '-x', 'DisplayName',  'SPA')
% hold on
% semilogx(list_eps, mean(RAE_fw, 2), '-s', 'DisplayName',  'FW')
% semilogx(list_eps, mean(RAE_fg, 2), '-o', 'DisplayName',  'FG')
% ylabel('Relative approximation error')
% legend
%
%
% figure();
% semilogx(list_eps, mean(MRSA_spa, 2), '-x', 'DisplayName',  'SPA')
% hold on
% semilogx(list_eps, mean(MRSA_fw, 2), '-s', 'DisplayName',  'FW')
% semilogx(list_eps, mean(MRSA_fg, 2), '-o', 'DisplayName',  'FG')
% ylabel('MRSA')
% legend


SNR = list_SNR;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultLineMarkerSize', 9)
set(groot, 'DefaultLineLineWidth', 1.4)
set(groot, 'DefaultAxesFontSize', 14);

figure();
plot(SNR, mean(RAE_fw, 2), '-s', 'DisplayName',  '\texttt{MERIT}')
hold on
plot(SNR, mean(RAE_fg, 2), '-o', 'DisplayName',  '\texttt{FastGradient}')
plot(SNR, mean(RAE_spa, 2), '-x', 'DisplayName',  '\texttt{SPA}')
xlabel('SNR')
ylabel('\texttt{Relative approximation error}', 'Interpreter',  'latex')
legend('Location',  'NorthWest')
% saveas(gcf, './results/gillis_setting_SNR/RAE.eps', 'epsc')

tmp = max([max(MRSA_spa(:)) max(MRSA_fw(:)) max(MRSA_fg(:))]);
figure();
plot(SNR, mean(MRSA_fw, 2)/tmp*100, '-s', 'DisplayName',  '\texttt{MERIT}')
hold on
plot(SNR, mean(MRSA_fg, 2)/tmp*100, '-o', 'DisplayName',  '\texttt{FastGradient}')
plot(SNR, mean(MRSA_spa, 2)/tmp*100, '-x', 'DisplayName',  '\texttt{SPA}')
ylabel('\texttt{MRSA}', 'Interpreter',  'latex')
xlabel('SNR')
legend
% saveas(gcf, './results/gillis_setting_SNR/MRSA.eps', 'epsc')

figure();
plot(SNR, mean(succ2_fw, 2), '-s', 'DisplayName',  '\texttt{MERIT}')
hold on
plot(SNR, mean(succ2_fg, 2), '-o', 'DisplayName',  '\texttt{FastGradient}')
plot(SNR, mean(succ2_spa, 2), '-x', 'DisplayName',  '\texttt{SPA}')
ylabel('\texttt{Overlapping rate}', 'Interpreter',  'latex')
xlabel('SNR')
legend('Location',  'NorthEast')
% saveas(gcf, './results/gillis_setting_SNR/overlapping.eps', 'epsc')

figure();
plot(SNR, mean(succ_fw, 2), '-s', 'DisplayName',  '\texttt{MERIT}')
hold on
plot(SNR, mean(succ_fg, 2), '-o', 'DisplayName',  '\texttt{FastGradient}')
plot(SNR, mean(succ_spa, 2), '-x', 'DisplayName',  '\texttt{SPA}')
ylabel('\texttt{success rate}', 'Interpreter',  'latex')
xlabel('SNR')
legend('Location',  'NorthEast')
% saveas(gcf, './results/gillis_setting_SNR/success_rate.eps', 'epsc')

