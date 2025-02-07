%% Clear all things
clc; clear; close all; path(pathdef);
addpath('~/code/matlab/common')
addpath('~/code/matlab/common/visualization/lines')
addpath('~/code/matlab/common/PGD')
addpath('~/code/matlab/common/prox_ops')

addpath('./synthetic_exp')
addpath('./topic_modeling/baselines/PrecondSPA')
addpath('./utils')

addpath('./fw_core/')

M = 50; L = 100; N = 10;

[X, pure_pixel_set, S] = generate_data(1, 3, M, N, L, 'SNR', 15);
optimalValue = N;

save('debug.mat', 'X', 'pure_pixel_set');
load('debug.mat', 'X', 'pure_pixel_set')

options = struct;
options.maxIters = 50;
options.lambda = 20;
options.objTol = -1;
options.N = N;
options.checkpointIterval = 1;
options.firstIterMethod = 'NONE';
options.verbose = 1;
options.innerVerboseInterval = 6000;
options.numThreads = 4;

fprintf('From MEX dense\n');
options.backend = 'matlab';
options.debug = true;

hparams = [1e-1 1e-2 1e-3];
objSequence = {};
for i=1:numel(hparams)
    options.epsilon = hparams(i);
    [C_hat2, Tracking] = fw(X, options);
    objSequence{i} = Tracking.InnerTracking.obj;
end

lines = objSequence;
output_dir = './results/mu_affect/';
names = {'1e-1', '1e-2', '1e-3'};
linewidth = 1.4;
markersize = 9;
linestyle_pools = {'-', '--', '-.'};
for i=1:numel(lines)
    plot(1:options.maxIters, lines{i}, sprintf('%s', linestyle_pools{i}), ...
        'DisplayName', names{i}, ...
        'LineWidth', linewidth,  ...
        'MarkerSize', markersize ...
        );
    hold on
end
xlabel('\mu');
ylabel('Objective function');
legend('Location', 'northeast');

set(gca, 'FontSize', 15);
    path_to_file = sprintf('%s/succ_rate.eps', output_dir);
    saveas(gcf, path_to_file, 'epsc')
    fprintf(sprintf('Exported to %s\n', path_to_file));

% plot_success_rate(1:100, objSequence, ...
%     names, './results/mu_affect/', '\mu', ...
%      'ylabel', 'Objective value');
%

