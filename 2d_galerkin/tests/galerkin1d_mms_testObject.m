% galerkin1d_mms_test.m
% Method of Manufactured Solutions (MMS) Test for GalerkinSolver1d Class
% Verifies the accuracy and convergence of the FEM solver by comparing
% numerical solutions against a known analytical solution.

clear; clc; close all;

%% 1. Define the True Solution and Compute Source Term

% Import Symbolic Toolbox
syms u_true_sym(x)

% (1) Specify the True Solution as a Symbolic Function
% Example: u_true(x) = sin(pi * x)
u_true_sym(x) = sin(pi * x);  % User can modify this as desired

% Physical Parameters
k = 1;  % Diffusion coefficient (k > 0)
r = 1; % Reaction coefficient (r > 0)

% (2) Compute the Source Term f(x) = -k * d²u_true/dx² + r * u_true using Symbolic Differentiation
f_sym = -k * diff(u_true_sym, x, 2) + r * u_true_sym;

% (3) Convert u_true and f_sym to MATLAB Function Handles
uTrue_func = matlabFunction(u_true_sym, 'Vars', x);
f_func = matlabFunction(f_sym, 'Vars', x);

% Store in Structures as per Existing Code
cofs.k = k;                 % Diffusion coefficient
cofs.r = r;                 % Reaction coefficient
cofs.f = @(x) f_func(x);    % Source term f(x)

% Display the Source Function for Verification
fprintf('Source Function f(x) = %s\n', char(f_sym));

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

% Preallocate for Convergence Rates
MMS.convergence_rate = zeros(num_tests-1, 1);

%% 3. Execute FEM Solver for Each Domain and Compute Errors

for i = 1:num_tests
    % Current Number of Finite Elements
    M_current = M_values(i);
    
    % Define the Domain Structure
    domain.M = M_current;                                     % Number of finite elements
    domain.L = 1;                                              % Length of the domain
    domain.x = linspace(0, domain.L, domain.M + 1)';         % Grid points (column vector)
    domain.h = domain.x(2) - domain.x(1);                     % Element size
    domain.N = domain.M - 1;                                   % Number of Interior Nodes (excluding boundaries)
    
    % Store the mesh size
    MMS.h(i) = domain.h;
    
    %% Instantiate and Run the GalerkinSolver1d Class
    
    % Create an instance of GalerkinSolver1d
    fem_solver = GalerkinSolver1d(domain, cofs);
    
    % Retrieve the Numerical Solution
    u_num = fem_solver.solution;  % Full solution including boundary nodes
    MMS.u_numeric{i} = u_num;
    
    %% Evaluate the True Solution on the Current Grid
    
    u_true_evaluated = uTrue_func(domain.x);
    MMS.u_true{i} = u_true_evaluated;
    
    %% Compute the L2 Error
    
    % Calculate the Error Vector
    error_vector = u_num - u_true_evaluated;
    
    % Compute the L2 Norm of the Error
    L2_error = sqrt(sum(error_vector.^2) * domain.h);
    MMS.errors(i) = L2_error;
    
    %% Display Progress
    fprintf('Test %d: M = %d, h = %.5f, L2 Error = %.5e\n', ...
        i, M_current, domain.h, L2_error);
end

%% 4. Compute Convergence Rates

for i = 1:num_tests-1
    rate = log(MMS.errors(i) / MMS.errors(i+1)) / log(MMS.h(i) / MMS.h(i+1));
    MMS.convergence_rate(i) = rate;
end

%% 5. Display Convergence Rates

%{
fprintf('\nConvergence Rates:\n');
for i = 1:num_tests-1
    fprintf('From M = %d to M = %d: Rate ≈ %.2f\n', ...
        MMS.Ms(i), MMS.Ms(i+1), MMS.convergence_rate(i));
end
%}

%% 6. Plotting the Errors vs. Mesh Size

%{
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
%}

%% 7. Summary of Results

% Create a Combined Table for Easy Viewing
Convergence_Rate = [NaN; MMS.convergence_rate];

% Construct the Combined Convergence Table
Convergence_Table = table(MMS.Ms', MMS.h, MMS.errors, Convergence_Rate, ...
    'VariableNames', {'M', 'h', 'L2_Error', 'Convergence_Rate'});

% Display the Combined Convergence Table
disp(' ');
disp('--- MMS Test Results ---');
disp(Convergence_Table);


%% 8. (Optional) Plot Numerical vs. True Solutions for Each Test

% Uncomment the following block to generate individual plots for each test

% for i = 1:num_tests
%     figure;
%     plot(MMS.u_numeric{i}, MMS.domain.x, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
%     hold on;
%     plot(MMS.u_true{i}, MMS.domain.x, 'r--', 'LineWidth', 2);
%     xlabel('Solution u(x)');
%     ylabel('Spatial Coordinate x');
%     title(sprintf('FEM Solution vs. True Solution (M = %d)', MMS.Ms(i)));
%     legend('Numerical FEM Solution', 'True Solution');
%     grid on;
%     hold off;
% end