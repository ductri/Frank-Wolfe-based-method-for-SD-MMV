%% Clear all things
clc; clear; close all; path(pathdef);
addpath('~/code/matlab/common')
addpath('~/code/matlab/common/visualization/lines')
addpath('~/code/matlab/common/PGD/')
addpath('~/code/matlab/common/prox_ops')

addpath('baselines/FGNSR/matlab')
addpath('baselines/PrecondSPA/')
addpath('fw_core')
addpath('utils')

%% UNIFORM A AND DIRICHLET S AND GAUSSIAN NOISE
M = 50; N = 100; K = 10;

global_tic = tic;
output_dir = './results/unknown_K/';
if ~isdir(output_dir)
    mkdir(output_dir)
    fprintf('Created output directory\n');
end
[X, pure_pixel_set, W, H, V] = generate_data(M, N, K, 'SNR', 20);
fw_tic = tic;
options = struct;
options.maxIters = 100;
options.verbose = 0;
options.lambda = 0.5;
options.debug = 0;
options.N = -1;
options.backend = 'mex';
options.dualGapThreshold=0;
options.epsilon = 1e-5;
[C_hat, fw_tracking] = fw(X, options);
[~, Lambda_hat] = maxk(vecnorm(C_hat, Inf, 2), K, 1);

figure();
linewidth = 1.4;
markersize = 6;
marker_pools = {'s', 'o', '^', 'x', 's', 'd', '^', 'v'};
linestyle_pools = {'-', '--', '-.'};

stem(vecnorm(C_hat, Inf, 2), 'bo', ...
        'LineWidth', linewidth,  ...
        'MarkerSize', markersize ...
        );
hold on
xline(pure_pixel_set, '--r');


% plot([1:N], vecnorm(C_hat, Inf, 2), 'bo', ...
%         'LineWidth', linewidth,  ...
%         'MarkerSize', markersize ...
%         );
% hold on
% xline(pure_pixel_set, '--b');

xlabel('index');
ylabel('||C(i, :)||_{\infty}');
ylim([0 1.0]);
set(gca, 'FontSize', 15);
saveas(gcf, sprintf('./%s/unknown_K.eps', output_dir), 'epsc')


