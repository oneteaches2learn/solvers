clear all; x = sym('x',[1 2],'real'); syms t; syms u; syms v;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% functions to project
u_fixed = exp(x(1) + x(2));
v_fixed = 10000;

% domain bounds
xLim_dom = [0 1];
yLim_dom = [0 1];

% mms parameters
p = 2;
base = 6;

% specify coefficients
k = 1;
lambda = 1;
r = lambda * u;
dr = lambda;

% specify BCs
bcTypes = 'NNNN';

% specify Neumann conditions
u_N = -(u^3 - v^3);
u_N = -(u^3 - v^3);
%u_N = v^4;
%u_N = u;
du_N = diff(u_N,u);


% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
auxfun = struct();

exp = [-5:1:5];
val = 10.^exp;
val = [-fliplr(val) 0 val];

%val = [-1.1:.1:0]

for v_fixed = val

    fprintf('2d Projection, v = %d\n',v_fixed)

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

    % solve
    fprintf(' Solving:'), tic
        prob = prob.solve;
    executionTime = toc; 
    fprintf(' %f s\n',executionTime)

    % print useful data
    fprintf('   avg(b_alpha) = %f\n', mean(full(prob.vectors.b_alpha_u)))
    fprintf('   avg(b_vol)   = %f\n', mean(full(prob.vectors.b_vol)))
    fprintf('   avg(M_r)     = %f\n', mean(full(prob.vectors.M_r)))
    fprintf('   avg(M_dr)    = %f\n', mean(full(prob.tensors.M_dr(:))))
    % plot
    fprintf(' Plotting\n\n')
    prob.plot;
    pause()

end


%{
% stash projection extras where the class expects them
prob.options.projection.u_fixed  = auxfun.u_fixed;
prob.options.projection.v_fixed  = auxfun.v_fixed;
prob.options.projection.alpha    = auxfun.alpha;
%}
