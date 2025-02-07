%% Clear all things
clc; clear; close all; path(pathdef);
addpath('baselines/PrecondSPA/')
addpath('~/code/matlab/common/prob_tools')
addpath('~/code/matlab/common/')
addpath('utils')
addpath('~/code/matlab/common/PGD/')
addpath('~/code/matlab/common/prox_ops')
addpath('fw_core')
addpath('baselines/FGNSR/matlab')
set(groot, 'DefaultLineMarkerSize', 9)
set(groot, 'DefaultLineLineWidth', 1.4)

M=50; N=200; K=40;
W = rand(M, K);
% W = dirichlet_rnd(ones(M, 1), K);
H = zeros(K, N);
H(:, N-K+1:end) = eye(K);
H(:, 1:N-K) = dirichlet_rnd(1*ones(K, 1), N-K);

pure_pixel_set = N-K+1:N;

spa_success_rate = [];
spa_success_count = [];
spa_mse = [];
spa_mse_l2 = [];
spa_Wmse_l2 = [];
spa_Wmse_l2_norm = [];

fw_success_rate = [];
fw_success_count = [];
fw_mse = [];
fw_mse_l2 = [];
fw_Wmse_l2 = [];
fw_Wmse_l2_norm = [];

fg_success_rate = [];
fg_success_count = [];
fg_mse = [];
fg_mse_l2 = [];
fg_Wmse_l2 = [];
fg_Wmse_l2_norm = [];

list_SNR = [2:1:16];
for j=1:numel(list_SNR)
    j
    for i=1:50
        Y = W*H;
        a = Y(:, 1:K);
        b = [Y(:, 1:K/2) Y(:, K+1:K+K/2)];
        norm(a-b, 'fro')^2/norm(W, 'fro')^2;
        e = mse_measure(a, b)/norm(W, 'fro')^2

        SNR = 10^(list_SNR(j)/10);
        noise = randn(size(Y)); 
        sigma2 = norm(Y, 'fro')^2 / M / N / SNR;
        noise = sqrt(sigma2)*noise;
        X = Y + noise;

        % ------------- SPA
        Lambda_hat = SPAselect(X, K);
        spa_success_count(j, i) = all(sort(Lambda_hat(:)) == sort(pure_pixel_set(:)));
        spa_success_rate(j, i) = numel(intersect(Lambda_hat, pure_pixel_set));
        W_hat = X(:, Lambda_hat);
        H_hat = computeApproxError(X, Lambda_hat);
        spa_mse(j, i) = sum(abs(X - W_hat*H_hat), 'all')/sum(abs(X), 'all');
        spa_mse_l2(j, i) = norm(X - W_hat*H_hat, 'fro')^2/norm(X, 'fro')^2;
        spa_Wmse_l2(j, i) = mse_measure(W_hat, W)/norm(W, 'fro')^2;
        spa_Wmse_l2_norm(j, i) = mse_measure_norm(W_hat, W);


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

        W_hat = X(:, Lambda_hat);
        H_hat = computeApproxError(X, Lambda_hat);
        fw_mse(j, i) = sum(abs(X - W_hat*H_hat), 'all')/sum(abs(X), 'all');
        fw_mse_l2(j, i) = norm(X - W_hat*H_hat, 'fro')^2/norm(X, 'fro')^2;
        fw_success_count(j, i) = all(sort(Lambda_hat(:)) == sort(pure_pixel_set(:)));
        fw_success_rate(j, i) = numel(intersect(Lambda_hat, pure_pixel_set));
        fw_Wmse_l2(j, i) = mse_measure(W_hat, W)/norm(W, 'fro')^2;
        fw_Wmse_l2_norm(j, i) = mse_measure_norm(W_hat, W);


        [C_hat_fg, K_fg, ~, ~, ~, ~, t] = fgnsr(X, K, 'maxiter', 150, ...
            'debug', 0, 'verbose', false);
        W_hat = X(:, K_fg);
        H_hat = computeApproxError(X, K_fg);
        fg_mse(j, i) = sum(abs(X - W_hat*H_hat), 'all')/sum(abs(X), 'all');
        fg_mse_l2(j, i) = norm(X - W_hat*H_hat, 'fro')^2/norm(X, 'fro')^2;
        fg_success_count(j, i) = all(sort(K_fg(:)) == sort(pure_pixel_set(:)));
        fg_success_rate(j, i) = numel(intersect(K_fg, pure_pixel_set));
        fg_Wmse_l2(j, i) = mse_measure(W_hat, W)/norm(W, 'fro')^2;
        fg_Wmse_l2_norm(j, i) = mse_measure_norm(W_hat, W);
    end
end
save('results/1_.mat')

figure();
plot(list_SNR, mean(spa_success_count, 2),'-x', 'DisplayName',  'SPA')
hold on
plot(list_SNR, mean(fw_success_count, 2),'-o', 'DisplayName',  'FW')
plot(list_SNR, mean(fg_success_count, 2), '-s', 'DisplayName',  'FG')
ylabel('success rate')
xlabel('SNR')
legend
saveas(gcf, 'results/1_success_rate.eps', 'epsc')

figure();
plot(list_SNR, mean(spa_success_rate, 2)/K,'-x', 'DisplayName',  'SPA')
hold on
plot(list_SNR, mean(fw_success_rate, 2)/K,'-o', 'DisplayName',  'FW')
plot(list_SNR, mean(fg_success_rate, 2)/K, '-s', 'DisplayName',  'FG')
ylabel('success rate 2')
xlabel('SNR')
legend
saveas(gcf, './results/1_success_rate2.eps', 'epsc')


figure();
plot(list_SNR, mean(spa_mse_l2, 2),'-x', 'DisplayName',  'SPA')
hold on
plot(list_SNR, mean(fw_mse_l2, 2),'-o', 'DisplayName',  'FW')
plot(list_SNR, mean(fg_mse_l2, 2), '-s', 'DisplayName',  'FG')
ylabel('NMSE')
xlabel('SNR')
legend
saveas(gcf, './results/1_mse_l2.eps', 'epsc')

figure();
plot(list_SNR, mean(spa_Wmse_l2, 2),'-x', 'DisplayName',  'SPA')
hold on
plot(list_SNR, mean(fw_Wmse_l2, 2),'-o', 'DisplayName',  'FW')
plot(list_SNR, mean(fg_Wmse_l2, 2), '-s', 'DisplayName',  'FG')
ylabel('W-NMSE')
xlabel('SNR')
legend
saveas(gcf, './results/1_Wmse_l2.eps', 'epsc')

figure();
plot(list_SNR, mean(spa_Wmse_l2_norm, 2),'-x', 'DisplayName',  'SPA')
hold on
plot(list_SNR, mean(fw_Wmse_l2_norm, 2),'-o', 'DisplayName',  'FW')
plot(list_SNR, mean(fg_Wmse_l2_norm, 2), '-s', 'DisplayName',  'FG')
ylabel('W-NMSE')
xlabel('SNR')
legend
saveas(gcf, './results/1_Wmse_l2_norm.eps', 'epsc')

