%% Clear all things
clc; clear; close all; path(pathdef);
addpath('~/code/matlab/common')
addpath('~/code/matlab/common/PGD')
addpath('~/code/matlab/common/prox_ops')

A = rand(10, 5);
x = dirichletRand(ones(1, 5), 1);
b = A*x;

x_hat = simplex_LS(A, b);
norm(x-x_hat)

