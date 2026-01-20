clear all; x = sym('x',[1 2],'real'); syms t; syms u; syms v;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% functions to project
u_fixed = exp(x(1) + x(2));
v_fixed = 2;

% domain bounds
xLim_dom = [0 1];
yLim_dom = [0 1];

% mms parameters
p = 2;
base = 4;

% specify coefficients
k = 1;
lambda = 1;
r = lambda * u;
dr = lambda;

% specify BCs
bcTypes = 'NNNN';

% specify Neumann conditions
u_N = - v * u^3;
du_N = diff(u_N,u);


% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
auxfun = struct();


% assemble domain
fprintf(' Contructing Domain:'), tic
    BCs = {u_N,u_N,u_N,u_N};
	dom = Domain2d(xLim_dom,yLim_dom);
    dom = GalerkinAssembler2d_frostbite.assembleBoundary(dom,bcTypes,BCs); 
    dom = GalerkinAssembler2d_frostbite.assembleNonlinearBoundary(dom,bcTypes,BCs);
    dom = dom.setMesh(p,base);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% assemble coefficients
fprintf(' Assembling Coefficients:'), tic
    % coefficients used by base assembler
    auxfun.cofs.k = matlabFunction(symfun(k,x));
    auxfun.cofs.r = matlabFunction(symfun(r,x));
    auxfun.cofs.dr_du = matlabFunction(symfun(dr,x));
    auxfun.f = @(x,y,t,U,V) 0;

    % projection-specific extras
    auxfun.params.u_fixed = matlabFunction(symfun(u_fixed,x));
    auxfun.params.v_fixed = v_fixed;
    auxfun.params.lambda  = 1.0;
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% assemble solver
fprintf(' Assembling Solver:'), tic
    prob = NewtonGalerkinSolver2d_projection(dom, auxfun);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% stash projection extras where the class expects them
prob.options.projection.u_fixed  = auxfun.u_fixed;
prob.options.projection.v_fixed  = auxfun.v_fixed;
prob.options.projection.alpha    = auxfun.alpha;
