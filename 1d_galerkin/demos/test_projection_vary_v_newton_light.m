% test_projection_vary_v_newton_light.m
%
% Lightweight Newton solver for the 1D P1 projection test problem
%   Omega=(0,1), mesh nodes 0,1/2,1
%   u_fixed(x)=x^2
%   alpha(u,v) = v*u  (Robin-type boundary flux)
%
% For this alpha, the discrete system is linear in w:
%   A(v,lambda) w = b(v,lambda)
% Newton's method for J(w)=A*w - b converges in one step in exact arithmetic.
%
% Use this as a "known-good" Newton reference before moving to truly
% nonlinear alpha(u,v).

%clear; clc;

%% parameters you might want to change
lambda       = 0;      % try 0 first (gives w=[0;1/4;1] independent of v)
newton_tol   = 1e-12;
newton_maxit = 25;

%% v sweep (same structure you used)
exponents_pos = 10.^(-7:1:7);
exponents_neg = -fliplr(exponents_pos);
exponents     = [exponents_neg exponents_pos];

%% fixed FE tensors for this mesh
K = [ 2 -2  0;
     -2  4 -2;
      0 -2  2];

M = [1/6  1/12 0;
     1/12 1/3  1/12;
     0    1/12 1/6];

B = diag([1 0 1]);   % boundary inner product in 1D: f(0)g(0)+f(1)g(1)

%% fixed parts of the RHS from u_fixed(x)=x^2
rhs_stiff = [-1/2; -1; +3/2];
rhs_mass  = [ 1/96; 7/48; 17/96];

%% Newton driver
fprintf('lambda = %.6g\n', lambda);
fprintf('%12s  %10s  %14s  %14s\n', 'v', 'iters', '||J(w)||_2', '||w-u_exact||_2');

%for v_fixed = exponents
for v_fixed = 10

    % assemble A and b for this v
    A = K + lambda*M - v_fixed*B;
    b = rhs_stiff + lambda*rhs_mass - v_fixed*[0;0;1];
    full(b)

    % exact (closed form) solution from the derived algebraic formula
    u_exact = nan(3,1);
    denom = 24 * (-lambda^2 + 8*lambda*v_fixed - 48*lambda + 96*v_fixed);
    if abs(denom) > 1e-14
        u_exact(1) =  lambda*(lambda + 48);
        u_exact(2) = -5*lambda^2 + 36*lambda*v_fixed - 240*lambda + 576*v_fixed;
        u_exact(3) = -23*lambda^2 + 192*lambda*v_fixed - 1104*lambda + 2304*v_fixed;
        u_exact = u_exact / denom;
    end

    % residual and Jacobian handles (Newton form)
    J  = @(w) A*w - b;
    DJ = @(w) A;  % constant for this linear alpha(u,v)

    % initial guess
    w0 = zeros(3,1);

    % run Newton
    [w,it,resnorm] = newton_solve(J, DJ, w0, newton_tol, newton_maxit);

    % error to exact (if available)
    if any(isnan(u_exact))
        err = nan;
    else
        err = norm(w - u_exact);
    end

    fprintf('%12.4e  %10d  %14.6e  %14.6e\n', v_fixed, it, resnorm, err);
end

%% ---------- local function: Newton solver ----------
function [w,it,resnorm] = newton_solve(J, DJ, w0, tol, maxit)
%NEWTON_SOLVE  Basic Newton iteration for solving J(w)=0.
%
%   [w,it,resnorm] = newton_solve(J, DJ, w0, tol, maxit)
%
% Inputs
%   J      residual handle: J(w)
%   DJ     Jacobian handle: DJ(w) (matrix)
%   w0     initial guess
%   tol    stopping tolerance on 2-norm of residual
%   maxit  maximum iterations
%
% Outputs
%   w       final iterate
%   it      number of iterations taken
%   resnorm final residual norm

w = w0;
res = J(w);
resnorm = norm(res);

it = 0;
while (resnorm > tol) && (it < maxit)
    A = DJ(w);

    % Solve A * delta = -J(w)
    delta = -(A \ res);

    % Update
    w = w + delta;

    it = it + 1;

    % Recompute residual
    res = J(w);
    resnorm = norm(res);
end

end
