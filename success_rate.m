%% Clear all things
clc; clear; close all; path(pathdef);
addpath('~/code/matlab/common')
addpath('~/code/matlab/common/visualization/lines')
addpath('~/code/matlab/common/PGD/')
addpath('~/code/matlab/common/prox_ops')
addpath('~/code/matlab/common/prob_tools')

addpath('baselines/FGNSR/matlab')
addpath('baselines/PrecondSPA/')
addpath('fw_core')
addpath('utils')

M=50; N=200; K=20;
hparams = [5:1:15];
num_trials = 50;
alpha = 0.5*ones(K, 1);

success_count = [];

Tracking = struct;
Tracking.lambda = -3;
Tracking.M = M;
Tracking.N = N;
Tracking.K = K;
Tracking.num_trials = num_trials;
Tracking.hparams = hparams;

global_tic = tic;
output_dir = sprintf( './results/SNR/%d-%d/', M, K);
if ~isdir(output_dir)
    mkdir(output_dir)
    fprintf('Created output directory\n');
end

fprintf('Running X2\n')
num_exps = numel(hparams);
for exp_ind=1:num_exps
    fprintf('Exp %d/%d -  Hparam: %f\n', exp_ind, num_exps, ...
        hparams(exp_ind))
    Tracking.fw(exp_ind).success_count = [];
    Tracking.fw0(exp_ind).success_count = [];
    Tracking.fg(exp_ind).success_count = [];
    Tracking.spa(exp_ind).success_count = [];

    for trial_ind=1:num_trials
        if mod(trial_ind, 10) == 0
            fprintf('Trial %d\n', trial_ind);
            fprintf('Total elapsed time: %.2f s\n', duration);
        end
            
        SNR = hparams(exp_ind);

        W = rand(M, K);
        H = zeros(K, N);
        H(:, 1:K) = eye(K);
        H(:, K+1:end) = dirichlet_rnd(alpha, N-K);
        Y = W*H;
        SNR = 10^(SNR/10);
        noise = randn(size(Y)); 
        sigma2 = sum(vecnorm(Y, 2, 1).^2) / M / N / SNR;
        noise = sqrt(sigma2)*noise;
        X = Y + noise;
        indices = randperm(N);
        X = X(:, indices);
        H = H(:, indices);
        r_pure_pixel_set = [];
        pure_pixel_set = 1:K;
        for ii=1:numel(pure_pixel_set)
            r_pure_pixel_set(end+1) = find(indices == pure_pixel_set(ii));
        end
        pure_pixel_set = r_pure_pixel_set;



        % ==============================
        % FRANK-WOLFE
        fw_tic = tic;
        options = struct;
        options.maxIters = 100;
        options.verbose = 0;
        options.lambda = Tracking.lambda;
        options.debug = 0;
        options.N = K;
        options.backend = 'mex';
        options.dualGapThreshold=0;
        options.epsilon = 1e-5;
        [C_hat, fw_tracking] = fw(X, options);
        [~, Lambda_hat] = maxk(vecnorm(C_hat, Inf, 2), K, 1);

        Tracking.fw(exp_ind).success_count(trial_ind) = ...
            all(sort(Lambda_hat(:)) == sort(pure_pixel_set(:)));
        Tracking.fw(exp_ind).duration(trial_ind) = toc(fw_tic);
        Tracking.fw(exp_ind).lambda(trial_ind) = fw_tracking.lambda;


        % ==============================
        options.lambda = 0;
        [C_hat, fw_tracking] = fw(X, options);
        [~, Lambda_hat] = maxk(vecnorm(C_hat, Inf, 2), K, 1);
        Tracking.fw0(exp_ind).success_count(trial_ind) = ...
            all(sort(Lambda_hat(:)) == sort(pure_pixel_set(:)));
        Tracking.fw0(exp_ind).duration(trial_ind) = toc(fw_tic);
        Tracking.fw0(exp_ind).lambda(trial_ind) = fw_tracking.lambda;



        % ==============================
        % FAST GRADIENT
        fg_tic = tic;
        [C_hat_fg, K_fg, ~, ~, ~, ~, t] = fgnsr(X, K, 'maxiter', 150, ...
            'debug', 0, 'verbose', false);
        Tracking.fg(exp_ind).success_count(trial_ind) = ...
            all(sort(K_fg(:)) == sort(pure_pixel_set(:)));
        Tracking.fg(exp_ind).duration(trial_ind) = toc(fg_tic);


        % ==============================
        % SPA
        i_tic = tic;
        Lambda_hat = SPAselect(X, K);
        Tracking.spa(exp_ind).success_count(trial_ind) = ...
            all(sort(Lambda_hat(:)) == sort(pure_pixel_set(:)));
        Tracking.spa(exp_ind).duration(trial_ind) = toc(i_tic);

        duration = toc(global_tic);
        Tracking.duration = duration;


        save(sprintf('%s/success_rate_%.2f.mat', output_dir, alpha(1)), 'Tracking')
    end
end

load(sprintf('%s/success_rate_%.2f.mat', output_dir, alpha(1)))
fprintf('Running X2\n')
fw = arrayfun(@(x) mean(x.success_count), Tracking.fw);
fg = arrayfun(@(x) mean(x.success_count), Tracking.fg);
spa = arrayfun(@(x) mean(x.success_count), Tracking.spa);
fw0 = arrayfun(@(x) mean(x.success_count), Tracking.fw0);

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultLineMarkerSize', 9)
set(groot, 'DefaultLineLineWidth', 1.4)
set(groot,'defaultAxesFontSize',20)

fw_color = [0 0.4470 0.7410];
fg_color = '#EDB120';
spa_color = [0.4940, 0.1840, 0.5560];

figure('DefaultAxesFontSize', 18);
plot(hparams, fw, '-s',  'DisplayName', '\texttt{MERIT}', 'Color', fw_color)
hold on
% plot(hparams, fw0, '--s',  'DisplayName', '\texttt{MERIT}(0)', 'Color', fw_color)
plot(hparams, fg, '-o',  'DisplayName', '\texttt{FastGradient}', 'Color', fg_color)
plot(hparams, spa, '-x',  'DisplayName', '\texttt{SPA}', 'Color', spa_color)

xlabel('SNR');
ylabel('success rate');
legend('Location', 'northwest', 'Interpreter', 'latex');
axis tight;
% set(gca, 'FontSize', 15);
saveas(gcf, sprintf('%s/harder-%d-%d-%.2f.eps', output_dir, M, K, alpha(1)), 'epsc')


