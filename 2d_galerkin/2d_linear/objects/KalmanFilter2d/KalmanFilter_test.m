runExperiment = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%actFunc = StolwijkActivation();
%sigma = actFunc.value;

%cofs = Coefficients;


    % DOMAIN INFORMATION
    %hand_small_fine
    %tyler_closed_medium
    %tyler_closed_fine
    tyler_closed_lc2

    % mirror the mesh across horizontal midline (if using tyler_closed_lc2)
    p = msh.POS;
    t = msh.TRIANGLES;
    %p(:,1) = -p(:,1);
    y0 = (min(p(:,1)) + max(p(:,1))) / 2;
    p(:,1) = 2*y0 - p(:,1);
    t(:,[2 3]) = t(:,[3 2]);
    msh.POS = p;
    msh.TRIANGLES = t;

    % rescale mesh
    msh.POS = msh.POS / max(max(msh.POS));
    msh.POS = msh.POS * 0.24;

    % Construct Domain Geometry
    dom_geo = Domain2d.domainFromGmsh(msh);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test of random variables
syms u v a b 

expr = a * (u - v) + b * u;

F = GaussianFunction(expr, [u, v]);

F = F.addParameter("a", GaussianParameter(10, 1));
F = F.addParameter("b", GaussianParameter(1000, 100));

[f1, params1] = F.sample()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test of random coefficients
c = @(x1,x2) 3.0e6 + 1.0e5 .* x1;

k = @(x1,x2) 0.45 + 0.05 .* x2;

f = 350;
uInit = 33;

auxRand = AuxFun_gaussian(c, k, f, uInit);

aux1 = auxRand.sample();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}


%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms x1 x2 a real

aRand = GaussianParameter(2.0, 0.25);

c = GaussianFunction(1 + a*cos(x1)^2, [x1, x2]);
c = c.addParameter("a", aRand);

k = 0.45;
f = 350;
uInit = 33;

auxRand = AuxFun_gaussian(c, k, f, uInit);

[aux1, samples1] = auxRand.sample();
[aux2, samples2] = auxRand.sample();

samples1.c
samples2.c

aux1.cofs.c(0.5, 0.7)
aux2.cofs.c(0.5, 0.7)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%value = f1(33, 37);

if runExperiment
    % rxndiffSolve_demo 
%    clearvars -except saveResults1 saveResults2 outputFolder plot_temps
    x = sym('x',[1 2],'real'); syms t; syms u; syms v;
    % USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % DOMAIN INFORMATION
    %hand_small_fine
    %tyler_closed_medium
    %tyler_closed_fine
    tyler_closed_lc2

    % mirror the mesh across horizontal midline (if using tyler_closed_lc2)
    p = msh.POS;
    t = msh.TRIANGLES;
    %p(:,1) = -p(:,1);
    y0 = (min(p(:,1)) + max(p(:,1))) / 2;
    p(:,1) = 2*y0 - p(:,1);
    t(:,[2 3]) = t(:,[3 2]);
    msh.POS = p;
    msh.TRIANGLES = t;

    % rescale mesh
    msh.POS = msh.POS / max(max(msh.POS));
    msh.POS = msh.POS * 0.24;

    % time stepping
    T  = 7200; 
    dt = 10; 

    % PDE INFORMATION
    % specify physical conditions
    u_body = 37;
    u_hand = 33;
    u_air = 0.5;

    % specify coefficients
    rho_blood = 1050;
    c_blood = 3617;
    c_vol = 3.5 * 10^3;
    c_hand = c_vol * 1000;
    c = c_hand;
    k = 0.42;
    r_const = rho_blood * c_blood;
    f = 0;

    % overwrite k and r_const with fake values to make me look good?
    %k = 2;
    %r_const = 0.5 * rho_blood * c_blood;

    % specify boundary conditions
    bcTypes_exterior = 'RNNN';
    bcTypes_interior = 'R';
    bc_air = 136 * (u - u_air);
    %bc_air = {136,u_air};
    bc_wrist = {10,v};
    BCs = {bc_wrist,bc_air,bc_air,bc_air,{0,0}};

    % specify initial condition
    u_o = u_hand;

    % ODE INFORMATION
    body_mass = 90;
    g = 350;
    %g = 800;
    v_o = u_body;
    c_body = c_vol;
    s_const = r_const * 0.01;
    constraints.vLower = NaN;
    constraints.vUpper = 37;
    order = 1;

    % ACTIVATION MODEL
    aType = 'Stolwijk';

    % BOUNDARY TYPE: 'nonlinear', 'linearized', or 'dirichlet'
    bcParm = 'linearized'; 


    % MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SETUP TRIAL
    fprintf('Solving Model\n')

    % build activation model
    if strcmp(aType,'C')
        r.r_const = r_const;
        r.r_activ = 1;
        r.r_activ_du = 0;
        kludge = GalerkinAssembler2d_frostbite.assembleCoefficients(c,k,f,u_o);
        kludge = GalerkinAssembler2d_frostbite.assembleReactionTerm(kludge,r.r_const,r.r_activ,r.r_activ_du);
        kludge_active = @(u,v) kludge.cofs.r_activ(0,0,0,u,v);
        data.cofs.s = @(u,v) s_const * kludge_active(u,v);
        s.s_func = data.cofs.s;
    else
        act.r_const = r_const;
        act.s_const = s_const;
        [r,s] = Coefficients.coupling_coefficients(aType,act);
    end

    % Construct Domain Geometry
    fprintf(' Contructing Domain Geometry:'), tic
        dom_geo = Domain2d.domainFromGmsh(msh);
    executionTime = toc; 
    fprintf(' %f s\n',executionTime)

    % Assemble PDE Data 
    fprintf(' Assembling PDE Data:'), tic
        bcTypes = [bcTypes_exterior, bcTypes_interior];
        dom = GalerkinAssembler2d_frostbite.assembleBoundary(dom_geo,bcTypes,BCs); 
        dom = GalerkinAssembler2d_frostbite.assembleNonlinearBoundary(dom,bcTypes,BCs);

        % replace boundary functions on edges 2-4
        bc = PhysiologicalBC;
        bc.uInf = u_air;
        for i = 2:4

            % use nonlinear BCs
            if strcmp(bcParm,'nonlinear')
                dom.boundary.edges(i).boundaryCondition = bc.BC;
                dom.boundary.edges(i).boundaryCondition_ddu = bc.BC_du;

            % use linearized BCs
            elseif strcmp(bcParm,'linearized')
                dom.boundary.edges(i).boundaryCondition = bc.BC_linearized;
                dom.boundary.edges(i).boundaryCondition_ddu = bc.BC_linearized_du;

            % use dirichlet BCs
            elseif strcmp(bcParm,'dirichlet')
                dom.boundary.edges(i).boundaryType = 'D';
                dom.boundary.edges(i).boundaryCondition = @(x1,x2,t,u,v) u_air;
                dom.boundary.edges(i).boundaryCondition_ddu = @(x1,x2,t,u,v) 0;
            end 
        end

        % continue processing domain
        dom.boundary = dom.boundary.setBoundaryNodeLists(dom.mesh);
        dom.boundary = dom.boundary.setFreeNodes(dom.mesh);
        dom = GalerkinAssembler2d_frostbite.assembleTimeStepping(dom,T,dt);
        auxfun = GalerkinAssembler2d_frostbite.assembleCoefficients(c,k,f,u_o);
        auxfun = GalerkinAssembler2d_frostbite.assembleReactionTerm(auxfun,r.r_const,r.r_activ,r.r_activ_du);
    executionTime = toc;
    fprintf(' %f s\n',executionTime)


    % Assemble ODE Data
    fprintf(' Assembling ODE Data:'), tic
        data.g = g;
        data.vInit = v_o;
        %r_activ = @(u,v) auxfun.cofs.r_activ(0,0,0,u,v);
        %data.cofs.s = @(u,v) s.s_const * s.s_activ(0,0,0,u,v);
        data.cofs.s = s.s_func;
        data.cofs.ds = @(u,v) s.s_const * s.s_activ_dv(0,0,0,u,v);
        data.cofs.c = c_body;
        data.time.T = T;
        data.time.dt = dt;
        data.constraints = constraints;
        options.order = order;
        ode = ODE(data,options);
    executionTime = toc;
    fprintf(' %f s\n',executionTime)

    % Solve
    fprintf(' Solving:'), tic
%        ode  = ODE(data,options);
%        kf   = KalmanFilter2d(dom,auxfun,ode);
    executionTime = toc; 
    fprintf(' %f s\n',executionTime)

end


%{
u = [25:.1:34];
%plot(u,sigma(u,37))
sigma_old = auxfun.cofs.r_activ;
omega_b = 1.1e-3;

hold on
plot(u,sigma_old(0,0,0,u,37))
plot(u,omega_b * sigma(u,37))
hold off
%}

%{
D = actFunc.D;
hold on
plot(u,cofs.DI(37-37,u-33))
plot(u,D(u,37))
hold off
%}

%{
C = actFunc.C;
hold on
plot(u,cofs.CS(37-37,u-33))
plot(u,C(u,37))
hold off
%}
