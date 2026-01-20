% initialize solution storage
ind = 1;
exponents_pos = 10.^(-7:1:7);
exponents_neg = -fliplr(exponents_pos);
exponents = [exponents_neg exponents_pos];

N = length(exponents);

for i = exponents

    % set u and v
    v_fixed = i;

    % compute exact solution
    u_exact_denom = 24 * (-lambda^2 + 8*lambda*v_fixed - 48*lambda + 96*v_fixed);
    u_exact(1,1) = lambda * (lambda + 48);
    u_exact(2,1) = -5 * lambda^2 + 36 * lambda * v_fixed - 240 * lambda + 576 * v_fixed;
    u_exact(3,1) = -23 * lambda^2 + 192 * lambda * v_fixed - 1104 * lambda + 2304 * v_fixed;
    u_exact = u_exact / u_exact_denom;

    % compute numerical solution
    K = [2 -2 0; -2 4 -2; 0 -2 2];
    M = [1/6 1/12 0; 1/12 1/3 1/12; 0 1/12 1/6];
    B = [1 0 0; 0 0 0; 0 0 1];
    b = [-1/2; -1; 3/2] + lambda * [1/96; 7/48; 17/96] - v_fixed * [0; 0; 1];

    % assemble
    A = K + lambda * M - v_fixed * B;

    % solve
    w = A \ b;

    % compute error
    %prob.solution' - u_exact
    u_exact - w
    pause()

end