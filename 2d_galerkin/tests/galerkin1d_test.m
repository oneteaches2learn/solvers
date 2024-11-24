% galerkin_1d.m
% Galerkin Finite Element Solver for the Poisson Equation:
%   -k * d^2u/dx^2 = f(x) on Omega = (0,1)
% with homogeneous Dirichlet boundary conditions:
%   u(0) = 0 and u(1) = 0

clear; clc; close all;

%% 1. Define Structures for Parameters and Domain

% Coefficients Structure
cofs.k = 1;                          % Diffusion coefficient (k > 0)
cofs.f = @(x) sin(pi * x);          % Source term f(x)

% Domain Structure
domain.M = 20;                       % Number of finite elements (total nodes = M + 1)
domain.L = 1;                        % Length of the domain
domain.x = linspace(0, domain.L, domain.M + 1)';  % Grid points (column vector)
domain.h = domain.x(2) - domain.x(1);             % Element size
domain.N = domain.M - 1;             % Number of Interior Nodes (excluding boundaries)

%% 2. Run Finite Element Method

% Call the FEM function with domain and coefficients
u = run_fem(domain, cofs);

%% 3. (Optional) Save Results

% Uncomment below lines to save the solution
% save('FEM_Solution.mat', 'u', 'domain', 'cofs');
