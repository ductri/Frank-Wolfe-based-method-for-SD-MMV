%% Clear all things
clc; clear; close all; path(pathdef);
addpath('baselines/PrecondSPA/')
addpath('~/code/matlab/common/prob_tools')
addpath('utils')
addpath('~/code/matlab/common/PGD/')
addpath('~/code/matlab/common/prox_ops')
addpath('fw_core')
addpath('baselines/FGNSR/matlab')

M=50; N=55; K=10;


RAE_fw = [];
RAE_fw0 = [];
RAE_spa = [];
RAE_fg = [];

succ_fw = [];
succ_fw0 = [];
succ_spa = [];
succ_fg = [];

succ2_fw = [];
succ2_fw0 = [];
succ2_spa = [];
succ2_fg = [];

MRSA_fw = [];
MRSA_fw0 = [];
MRSA_spa = [];
MRSA_fg = [];

list_SNR = [2:1:30];
% list_SNR = [2:3:30];

t0 = tic;
for exp_ind=1:numel(list_SNR)
    fprintf('%d/%d - Elapsed time: %f \n', exp_ind, numel(list_SNR), toc(t0));

    for trial_ind=1:50
        W = rand(M, K);
        W = W./sum(W, 1);

        H = zeros(K, N);
        H(:, 1:K) = eye(K);
        col = K+1;
        for i=1:K
            for ii=i+1:K
                H(i, col) = 0.5;
                H(ii, col) = 0.5;
                col = col+1;
            end
        end
        Y = W*H;
        pure_pixel_set = 1:K;

        w_bar = mean(W, 2);
        original_noise = zeros(M, N);
        original_noise(:, K+1:end) = Y(:, K+1:end) - w_bar;

        noise_energy = sqrt(norm(Y, 'fro')^2/10^(list_SNR(exp_ind)/10));
        noise = original_noise/norm(original_noise, 'fro')*noise_energy;
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

        % % ------------- FRANK-WOLFE tunning
        hyperparams_fw = linspace(0, fw_tracking.lambda*2, 10);
        for param_ind=1:numel(hyperparams_fw)
            options.lambda = hyperparams_fw(param_ind);
            [C_hat, fw_tracking] = fw(X, options);
            [~, Lambda_hat] = maxk(vecnorm(C_hat, Inf, 2), K, 1);
            W_hat = X(:, Lambda_hat);
            H_hat = computeApproxError(X, Lambda_hat);
            RAE_fw_tunning(exp_ind, trial_ind, param_ind) = 1-norm(X - W_hat*H_hat, 'fro')/norm(X, 'fro');
            MRSA_fw_tunning(exp_ind, trial_ind, param_ind) = compute_MRSA(W_hat, W);
            succ_fw_tunning(exp_ind, trial_ind, param_ind) = all(sort(Lambda_hat(:)) == sort(pure_pixel_set(:)));
            succ2_fw_tunning(exp_ind, trial_ind, param_ind) = numel(intersect(Lambda_hat, pure_pixel_set))/K;
        end


        options.lambda = 0;
        [C_hat, fw_tracking] = fw(X, options);
        [~, Lambda_hat] = maxk(vecnorm(C_hat, Inf, 2), K, 1);
        W_hat = X(:, Lambda_hat);
        H_hat = computeApproxError(X, Lambda_hat);
        RAE_fw0(exp_ind, trial_ind) = 1-norm(X - W_hat*H_hat, 'fro')/norm(X, 'fro');
        MRSA_fw0(exp_ind, trial_ind) = compute_MRSA(W_hat, W);
        succ_fw0(exp_ind, trial_ind) = all(sort(Lambda_hat(:)) == sort(pure_pixel_set(:)));
        succ2_fw0(exp_ind, trial_ind) = numel(intersect(Lambda_hat, pure_pixel_set))/K;


        [C_hat_fg, K_fg, ~, ~, ~, ~, t] = fgnsr(X, K, 'maxiter', 150, ...
            'debug', 0, 'verbose', false);
        succ_fg(exp_ind, trial_ind) = all(sort(K_fg(:)) == sort(pure_pixel_set(:)));
        succ2_fg(exp_ind, trial_ind) = numel(intersect(K_fg, pure_pixel_set))/K;

        W_hat = X(:, K_fg);
        H_hat = computeApproxError(X, K_fg);
        RAE_fg(exp_ind, trial_ind) = 1-norm(X - W_hat*H_hat, 'fro')/norm(X, 'fro');
        MRSA_fg(exp_ind, trial_ind) = compute_MRSA(W_hat, W);

        hyperparams_fg = linspace(0, t.mu*2, 10);
        for param_ind=1:numel(hyperparams_fg)
            [C_hat_fg, K_fg, ~, ~, ~, ~, t] = fgnsr(X, K, 'maxiter', 150, ...
                'debug', 0, 'verbose', false, 'mu', hyperparams_fg(param_ind));
            W_hat = X(:, K_fg);
            H_hat = computeApproxError(X, K_fg);
            RAE_fg_tunning(exp_ind, trial_ind, param_ind) = 1-norm(X - W_hat*H_hat, 'fro')/norm(X, 'fro');
            MRSA_fg_tunning(exp_ind, trial_ind, param_ind) = compute_MRSA(W_hat, W);
            succ_fg_tunning(exp_ind, trial_ind, param_ind) = all(sort(K_fg(:)) == sort(pure_pixel_set(:)));
            succ2_fg_tunning(exp_ind, trial_ind, param_ind) = numel(intersect(K_fg, pure_pixel_set))/K;
        end

    end
end

% save('./results/fw_vs_fg/gillis_setting.mat')
load('./results/fw_vs_fg/gillis_setting.mat')


SNR = list_SNR;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultLineMarkerSize', 9)
set(groot, 'DefaultLineLineWidth', 1.4)
% set(groot, 'DefaultAxesFontSize', 14);

fw_color = [0 0.4470 0.7410];
fg_color = '#EDB120';
spa_color = [0.4940, 0.1840, 0.5560];

figure('DefaultAxesFontSize', 18);
plot(SNR, mean(RAE_fw, 2), '-s', 'DisplayName',  '\texttt{MERIT}', 'Color', fw_color)
hold on
% plot(SNR, max(mean(RAE_fw_tunning, 2), [], 3), 'b--s', 'DisplayName',  '\texttt{MERIT-tuning}', 'Color', fw_color)
% plot(SNR, mean(RAE_fw0, 2), '--s',  'DisplayName',  '\texttt{MERIT}(0)', 'Color', fw_color);
plot(SNR, mean(RAE_fg, 2), '-o', 'DisplayName',  '\texttt{FastGradient}', 'Color', fg_color)
% plot(SNR, max(mean(RAE_fg_tunning, 2), [], 3), '--o', 'DisplayName',  '\texttt{FastGradient-tuning}' , 'Color', fg_color)
plot(SNR, mean(RAE_spa, 2), '-x', 'DisplayName',  '\texttt{SPA}', 'Color', spa_color)
xlabel('SNR')
ylabel('\texttt{RAE}', 'Interpreter',  'latex')
legend('Location',  'NorthWest', 'Interpreter',  'latex')
saveas(gcf, './results/fw_vs_fg/RAE_include_spa.eps', 'epsc')


figure('DefaultAxesFontSize', 16);
plot(SNR, mean(RAE_fw, 2), '-s', 'DisplayName',  '\texttt{MERIT}', 'Color',  fw_color)
hold on
plot(SNR, max(mean(RAE_fw_tunning, 2), [], 3), '--s', 'DisplayName',  '\texttt{MERIT-tuning}', 'Color',  fw_color)
plot(SNR, mean(RAE_fg, 2), '-o', 'DisplayName',  '\texttt{FastGradient}', 'Color', fg_color)
plot(SNR, max(mean(RAE_fg_tunning, 2), [], 3), '--o', 'DisplayName',  '\texttt{FastGradient-tuning}' , 'Color', fg_color)
% plot(SNR, mean(RAE_spa, 2), '-x', 'DisplayName',  '\texttt{SPA}')
xlabel('SNR')
ylabel('\texttt{RAE}', 'Interpreter',  'latex')
legend('Location',  'NorthWest', 'Interpreter',  'latex')
saveas(gcf, './results/fw_vs_fg/RAE.eps', 'epsc')



tmp = max([max(MRSA_spa(:)) max(MRSA_fw(:)) max(MRSA_fg(:))]);
figure('DefaultAxesFontSize', 18);
plot(SNR, mean(MRSA_fw, 2)/tmp*100, '-s', 'DisplayName',  '\texttt{MERIT}', 'Color', fw_color)
hold on
% plot(SNR, mean(MRSA_fw0, 2)/tmp*100, '--s', 'DisplayName',  '\texttt{MERIT}(0)', 'Color', fw_color)
% plot(SNR, min(mean(MRSA_fw_tunning, 2), [], 3)/tmp*100, 'b--s', 'DisplayName',  '\texttt{MERIT-tuning}', 'Color', fw_color)
plot(SNR, mean(MRSA_fg, 2)/tmp*100, '-o', 'DisplayName',  '\texttt{FastGradient}', 'Color', fg_color)
% plot(SNR, min(mean(MRSA_fg_tunning, 2), [], 3)/tmp*100, '--o', 'DisplayName',  '\texttt{FastGradient-tuning}', 'Color', fg_color)
plot(SNR, mean(MRSA_spa, 2)/tmp*100, '-x', 'DisplayName',  '\texttt{SPA}', 'Color', spa_color)
ylabel('\texttt{MRSA}', 'Interpreter',  'latex')
xlabel('SNR')
legend
saveas(gcf, './results/fw_vs_fg/MRSA_include_spa.eps', 'epsc')


tmp = max([max(MRSA_spa(:)) max(MRSA_fw(:)) max(MRSA_fg(:))]);
figure('DefaultAxesFontSize', 16);
plot(SNR, mean(MRSA_fw, 2)/tmp*100, '-s', 'DisplayName',  '\texttt{MERIT}', 'Color', fw_color)
hold on
plot(SNR, min(mean(MRSA_fw_tunning, 2), [], 3)/tmp*100, '--s', 'DisplayName',  '\texttt{MERIT-tuning}', 'Color', fw_color)
plot(SNR, mean(MRSA_fg, 2)/tmp*100, '-o', 'DisplayName',  '\texttt{FastGradient}', 'Color', fg_color)
plot(SNR, min(mean(MRSA_fg_tunning, 2), [], 3)/tmp*100, '--o', 'DisplayName',  '\texttt{FastGradient-tuning}', 'Color', fg_color)
% plot(SNR, mean(MRSA_spa, 2)/tmp*100, '-x', 'DisplayName',  '\texttt{SPA}')
ylabel('\texttt{MRSA}', 'Interpreter',  'latex')
xlabel('SNR')
legend
saveas(gcf, './results/fw_vs_fg/MRSA.eps', 'epsc')




% figure();
% plot(SNR, mean(succ2_spa, 2), '-x', 'DisplayName',  '\texttt{SPA}')
% hold on
% plot(SNR, mean(succ2_fw, 2), '-s', 'DisplayName',  '\texttt{MERIT}')
% plot(SNR, mean(succ2_fg, 2), '-o', 'DisplayName',  '\texttt{FastGradient}')
% ylabel('\texttt{Overlapping rate}', 'Interpreter',  'latex')
% xlabel('SNR')
% legend('Location',  'NorthWest')


figure('DefaultAxesFontSize', 18);
plot(SNR, mean(succ_fw, 2), '-s', 'DisplayName',  '\texttt{MERIT}', 'Color', fw_color)
hold on
% plot(SNR, mean(succ_fw0, 2), '--s', 'DisplayName',  '\texttt{MERIT}(0)', 'Color', fw_color)
% plot(SNR, max(mean(succ_fw_tunning, 2), [], 3), 'b--s', 'DisplayName',  '\texttt{MERIT-tuning}')
plot(SNR, mean(succ_fg, 2), '-o', 'DisplayName',  '\texttt{FastGradient}', 'Color', fg_color)
% plot(SNR, max(mean(succ_fg_tunning, 2), [], 3), '--o', 'DisplayName',  '\texttt{FastGradient-tuning}', 'Color', fg_color)
plot(SNR, mean(succ_spa, 2), '-x', 'DisplayName',  '\texttt{SPA}', 'Color', spa_color)
ylabel('\texttt{success rate}', 'Interpreter',  'latex')
xlabel('SNR')
legend('Location',  'SouthEast')
saveas(gcf, './results/fw_vs_fg/success_rate_include_spa.eps', 'epsc')



figure('DefaultAxesFontSize', 16);
plot(SNR, mean(succ_fw, 2), '-s', 'DisplayName',  '\texttt{MERIT}', 'Color', fw_color)
hold on
plot(SNR, max(mean(succ_fw_tunning, 2), [], 3), '--s', 'DisplayName',  '\texttt{MERIT-tuning}', 'Color', fw_color)
plot(SNR, mean(succ_fg, 2), '-o', 'DisplayName',  '\texttt{FastGradient}', 'Color', fg_color)
plot(SNR, max(mean(succ_fg_tunning, 2), [], 3), '--o', 'DisplayName',  '\texttt{FastGradient-tuning}', 'Color', fg_color)
% plot(SNR, mean(succ_spa, 2), '-x', 'DisplayName',  '\texttt{SPA}')
ylabel('\texttt{success rate}', 'Interpreter',  'latex')
xlabel('SNR')
legend('Location',  'SouthEast')
saveas(gcf, './results/fw_vs_fg/success_rate.eps', 'epsc')




