%% Clear all things
clc; clear; close all; path(pathdef);
addpath('~/code/matlab/common')
addpath('~/code/matlab/common/prox_ops')

addpath('./fw_core/')

M=100; N=500; K=40;
[X, pure_pixel_set, W, H, V] = generate_data(M, N, K, 'SNR', 50);
delta = max(vecnorm(V, 2, 1));
lambda = 1e-3;
Gamma = max(vecnorm(W, 2, 1));
Alpha = norm(W(:, 1) - W(:, 2:end)*ones(K-1, 1)/(K-1));
Beta = (sqrt(4*1*(1-K/N)*norm(V, 'fro')^2 + 2*lambda*K) + 2*delta)/Alpha / 0.5

options = struct;
options.maxIters = 30;
options.lambda = lambda;
options.objTol = -1;
options.N = N;
options.checkpointIterval = 1;
options.firstIterMethod = 'NONE';
options.verbose = 1;
options.innerVerboseInterval = 6000;
options.numThreads = 4;

fprintf('From MEX dense\n');
tic
options.backend = 'mex';
options.epsilon = 1e-6;
[C_hat, Tracking] = fw(X, options);
toc

fprintf('\n');
sort(pure_pixel_set)
[v, lambdaHat] = maxk(vecnorm(C_hat, Inf, 2), K, 1);
sort(lambdaHat)'

% bar(vecnorm(X - X*C_hat, 2, 1))
