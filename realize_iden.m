%% Clear all things
clc; clear; close all; path(pathdef);
addpath('~/code/matlab/common')
addpath('~/code/matlab/common/prox_ops')

addpath('./fw_core/')
addpath('./utils')
addpath('~/code/matlab/common/PGD')

M=100; N=200; K=5;

W = 1*rand(M, K);
% W = 1000*dirichlet_rnd(100*ones(1, M), K);
% W(W<0.7) = 0;

H = zeros(K, N);
H(:, 1:K) = eye(K);
H(:, K+1:end) = dirichlet_rnd(0.5*ones(1, K), N-K);
Y = W*H;
SNR = 40; SNR = 10^(SNR/10);
noise = randn(size(Y)); 
sigma2 = sum(vecnorm(Y, 2, 1).^2) / M / N / SNR;
V = sqrt(sigma2)*noise;
X = Y + V;

indices = randperm(N);
X = X(:, indices);
H = H(:, indices);
r_pure_pixel_set = [];
pure_pixel_set = 1:K;
for i=1:numel(pure_pixel_set)
    r_pure_pixel_set(end+1) = find(indices == pure_pixel_set(i));
end
pure_pixel_set = r_pure_pixel_set;



delta = max(vecnorm(V, 2, 1));
lambda = norm(V, 'fro')^2*10
Gamma = max(vecnorm(W, 2, 1));

Alpha = min([simplex_LS(W(:, 2:end), W(:, 1)), ...
    simplex_LS(W(:, [1, 3, 4, 5]), W(:, 2)), ...
    simplex_LS(W(:, [1, 2, 4, 5]), W(:, 3)), ...
    simplex_LS(W(:, [1, 2, 3, 5]), W(:, 4)), ...
    simplex_LS(W(:, [1, 2, 3, 4]), W(:, 5)), ...
    ])

Beta = (sqrt(4*1*(1-K/N)*norm(V, 'fro')^2 + 2*lambda*K) + 2*delta)/Alpha / 0.5
noise = norm(V, 'fro')^2

