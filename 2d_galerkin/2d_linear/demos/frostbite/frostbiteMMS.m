%NOTE: $r(u,v) = u*v*(u - v)$, $s(u,v) = v*u*(v - u)$

% EXPERIMENT 1
runTrials = 1;
runErrors = 1;

x = sym('x',[1 2],'real'); syms t; syms u; syms v;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOMAIN DATA
% domain bounds
xLim_dom = [0 1];
yLim_dom = [0 1];

% Y bounds
xLim_Y = [0 1];
yLim_Y = [0 1];

% mesh parameters
base = 2;

% TIME STEPPING DATA
T  = 1;
eq.atEq = "nobreak";
eq.tolerance = 1e-5;

% PDE INFORMATION
c = 1;
k = 1;
r_const = 1;
phi = sin(pi * x(1)) * sin(pi * x(2));
rho = (t + 1);
f = rho^2 * exp(t) * phi^2 + (1 + 2 * pi^2 * rho - rho * exp(2*t)) * phi;

% specify BCs
bTypes = 'DDDD';
bTypes2 = 'D';

% specify Dirichlet conditions
u_D = 0;

% specify Neumann conditions
u_N = 0;

% specify Robin conditions
alpha = 1;
u_R = 1;

% specify initial condition
u_o = sin(pi * x(1)) * sin(pi * x(2));

% ODE INFORMATION
c_body = 1;
s_const = 1;
%g = @(t)(2 * exp(t) - 4 * (t + 1) / pi^2);
g = @(t)(exp(t) .* (1 + exp(t) - (t + 1)*4/pi^2));
v_o = 1;

% ODE constraints
constraints.vLower = NaN;
constraints.vUpper = NaN;
order = 1;


% TRUE SOLUTIONS
uTrue = @(x1,x2,t)(sin(pi * x1) .* sin(pi * x2) .* (t + 1));
vTrue = @(t)(exp(t));

if runTrials
    % SETUP TRIAL
    fprintf('Solving Model 2\n')

    % ramp activation model
    r_activ = @(x1,x2,t,u,v)(u*v);
    r_activ_du = @(x1,x2,t,u,v)(v);

    % assemble domain
    fprintf(' Contructing Domain:'), tic
        %inc = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
        %inc = Inclusion2d_square(xLim_Y,yLim_Y,incRatio);
        %dom = GalerkinAssembler2d_heat.assembleDomainGeometry(xLim_dom,yLim_dom,inc,eps);
        dom_geo = GalerkinAssembler2d_rxndiff.assembleDomainGeometry(xLim_dom,yLim_dom);
        %dom = dom.add_yline;
        %dom = dom.inclusionsON;
        dom_geo = GalerkinAssembler2d_rxndiff.assembleBoundary(dom_geo,bTypes,u_D,u_N,alpha,u_R,bTypes2); 
    executionTime = toc; 
    fprintf(' %f s\n',executionTime)

    % Assemble PDE Data 
    fprintf(' Assembling PDE Data:'), tic
        auxfun = GalerkinAssembler2d_frostbite.assembleCoefficients(c,k,f,u_o);
        auxfun = GalerkinAssembler2d_frostbite.assembleReactionTerm(auxfun,r_const,r_activ,r_activ_du);
    executionTime = toc;
    fprintf(' %f s\n',executionTime)

    % Assemble ODE Data
    fprintf(' Assembling ODE Data:'), tic
        data.g = g;
        data.vInit = v_o;
        data.cofs.s = @(u,v) s_const * v;
        data.cofs.c = c_body;
        data.time.T = T;
        %data.time.dt = dt;
        data.constraints = constraints;
        options.order = order;
    executionTime = toc;
    fprintf(' %f s\n',executionTime)

    for i = 1:4

        % Update Mesh Data
        p_space = 3 + (i-1)
        %p_time  = 3 + 2*(i-1)
        p_time  = 3 + (i-1)
        p = p_space;
        dt = 2^-p_time;
        data.time.dt = dt;
        dom = GalerkinAssembler2d_rxndiff.assembleTimeStepping(dom_geo,T,dt,eq);
        dom = GalerkinAssembler2d_rxndiff.assembleMesh(dom,p,base); 

        % Solve
        fprintf(' Solving:'), tic
            ode = ODE(data,options);
            prob{i} = CoupledNewtonSolver2d_frostbite(dom,auxfun,ode);
        executionTime = toc; 
        fprintf(' %f s\n',executionTime)
    end
end

if runErrors

    % store variables
    trials = length(prob);

    % Store errors
    errors = zeros(1,trials);
    ratios = zeros(1,trials-1);
    orders = zeros(1,trials-1);

    % Compute error
    fprintf('Computing Errors:\n')
    %uTrue = self.auxFunctions.uTrue;

    for i = 1:trials

        fprintf(' Trial: '); tic;
        % compute quadrature (on each timestep, if time-varying)
        sol = prob{i}.solution;
        err_pde = (prob{i}.domain.L2err_threePointQuadrature_nodal(sol,uTrue)').^2;
        %err_pde = (prob{i}.domain.L2err_centroidQuadrature_nodal(sol,uTrue)').^2;
        %err_pde = (prob{i}.domain.L2err_nodalQuadrature(sol,uTrue)').^2;
        err_ode = abs(prob{i}.ODE.solution - vTrue(prob{i}.ODE.tGrid)).^2;
        err = sqrt(err_pde + err_ode);

        % optional: print largest error for each component
        max_pde = max(sqrt(err_pde));
        max_ode = max(sqrt(err_ode));
        
        % compute L_infty error across all time steps
        errors(i) = max(err);

        executionTime = toc;
        fprintf('%f s\n',executionTime) 
    end

    % Compute ratios and orders
    for i = 2:trials
        ratios(i-1) = errors(i-1)/errors(i);
        orders(i-1) = log(ratios(i-1)) / log(base);
    end

    orders

end