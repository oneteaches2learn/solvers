% galerkin1d_mms_test.m
% Method of Manufactured Solutions (MMS) Test for Galerkin Finite Element Solver
% Verifies the accuracy and convergence of the FEM solver by comparing
% numerical solutions against a known analytical solution.

clear; clc; close all;

%% 1. Define the True Solution and Compute Source Term

% Import Symbolic Toolbox
syms u_true(x)

% (1) Specify the True Solution as a Symbolic Function
% Example: u_true(x) = sin(pi * x)
u_true(x) = sin(pi * x);  % User can modify this as desired

% Physical Parameters
k = 1;  % Diffusion coefficient (k > 0)

% (2) Compute the Source Term f(x) = -k * d²u_true/dx² using Symbolic Differentiation
f_sym = -k * diff(u_true, x, 2);

% (3) Convert u_true and f_sym to MATLAB Function Handles
uTrue_func = matlabFunction(u_true, 'Vars', x);
f_func = matlabFunction(f_sym, 'Vars', x);

% Store in Structures as per Existing Code
cofs.k = k;             % Diffusion coefficient
cofs.f = @(x) f_func(x);  % Source term f(x)

%% 2. Define a Set of Domain Structures for MMS Test

% (4) Create Four Domain Structures with M = 10, 100, 1000, 10000
M_values = [10, 100, 1000, 10000];
num_tests = length(M_values);

% Preallocate Structures for Efficiency
MMS.Ms = M_values;
MMS.errors = zeros(num_tests, 1);
MMS.h = zeros(num_tests, 1);
MMS.u_numeric = cell(num_tests, 1);
MMS.u_true = cell(num_tests, 1);

% Loop Over Each Domain to Execute FEM and Compute Errors
for i = 1:num_tests
    % Current Number of Finite Elements
    M_current = M_values(i);
    
    % Define the Domain Structure
    domain.M = M_current;                       % Number of finite elements
    domain.L = 1;                                % Length of the domain
    domain.x = linspace(0, domain.L, domain.M + 1)';  % Grid points (column vector)
    domain.h = domain.x(2) - domain.x(1);               % Element size
    domain.N = domain.M - 1;                         % Number of Interior Nodes
    
    % Store the mesh size
    MMS.h(i) = domain.h;
    
    % (5) Execute the Finite Element Method
    u_num = run_fem(domain, cofs);  % Numerical Solution including boundary nodes
    MMS.u_numeric{i} = u_num;
    
    % Evaluate the True Solution on the Current Grid
    u_true_evaluated = uTrue_func(domain.x);
    MMS.u_true{i} = u_true_evaluated;
    
    % Compute the L2 Error (Assuming Uniform Mesh)
    % Exclude boundary nodes for comparison
    error_vector = u_num - u_true_evaluated;
    L2_error = sqrt(sum(error_vector.^2) * domain.h);
    MMS.errors(i) = L2_error;
    
    fprintf('Test %d: M = %d, h = %.5f, L2 Error = %.5e\n', i, M_current, domain.h, L2_error);
end

%% 3. Compute Convergence Rates

% Initialize Convergence Rates Vector
MMS.convergence_rate = zeros(num_tests-1, 1);

for i = 1:num_tests-1
    rate = log(MMS.errors(i)/MMS.errors(i+1)) / log(MMS.h(i)/MMS.h(i+1));
    MMS.convergence_rate(i) = rate;
end

% Display Convergence Rates
fprintf('\nConvergence Rates:\n');
for i = 1:num_tests-1
    fprintf('From M = %d to M = %d: Rate ≈ %.2f\n', ...
        MMS.Ms(i), MMS.Ms(i+1), MMS.convergence_rate(i));
end

%% 4. (Optional) Plotting the Errors vs. Mesh Size

figure;
loglog(MMS.h, MMS.errors, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
xlabel('Mesh Size h');
ylabel('L2 Error');
title('Convergence of FEM Solver using MMS');
grid on;

% Fit a line to the log-log data to estimate the convergence rate
p = polyfit(log(MMS.h), log(MMS.errors), 1);
hold on;
loglog(MMS.h, exp(p(2)) * MMS.h.^p(1), 'r--', 'LineWidth', 2);
legend('Numerical Errors', sprintf('Fit: Slope = %.2f', p(1)));
hold off;

%% 5. Summary of Results

% Create a Table for Easy Viewing
if num_tests >1
    Rates = [MMS.errors; NaN];
    Rates = reshape([MMS.errors(1:end-1), MMS.errors(2:end)], [],1);
    Rates = MMS.convergence_rate;
end

Convergence_Table = table(MMS.Ms', MMS.h', MMS.errors, ...
    'VariableNames', {'Number_of_Elements_M', 'Mesh_Size_h', 'L2_Error'});

if num_tests >1
    Rates_Table = table((1:num_tests-1)', MMS.convergence_rate, ...
        'VariableNames', {'Test_Pair', 'Convergence_Rate'});
end

disp('--- MMS Test Results ---');
disp(Convergence_Table);

if num_tests >1
    disp('--- Convergence Rates ---');
    disp(Rates_Table);
end