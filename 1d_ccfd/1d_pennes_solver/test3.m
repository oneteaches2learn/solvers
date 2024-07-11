clear variables

% DATA
% Coefficients and Source
perm_fun   = @(x)(0.5 + 0*x);
pennes_fun = @(x)(1 + 0*x);
pennes_cof = 1;
source_fun = @(x)(0 + 0*x);
u_cutoff   = 10;

% Time Data
tau = 0.01;
nt  = 200;
init_fun = @(x)(37 + 0*x);

% Domain Data
a = 0; b = 1; M = 21;
cellEdges = linspace(a,b,M);

% Boundary Data
u_a = 0;
u_b = 37;

% ODE data
alpha = 0.1;
mu = 0;
v_init = 37;
v_coma = 32;
v_star = 37;


% SOLVER
% Construct domain and time-stepping
dom  = Domain(cellEdges);
time = TimeStepping(tau,nt); 

% Construct coefficients
init   = InitialCondition(dom,init_fun);
perm   = Permeability(dom,perm_fun);
source = Source(dom,perm_fun);
pennes_cof = Coefficient(dom,pennes_fun);
trans  = Transmissibility(dom,perm);
laplacian = DiscreteLaplacian(trans);

% Construct Boundary Conditions
left     = DirichletBoundary(u_a,"L",trans);
right    = DirichletBoundary(u_b,"R",trans);
boundary = BoundaryConditions(left,right,dom);

M  = dom.cellNum;
dx = dom.cellWidths;
tau = time.stepSize;
nt  = time.stepNum;

% PENNES SOLVER
% Initialize data storage
W = zeros(dom.cellNum+1,time.stepNum);
U = zeros(dom.cellNum,time.stepNum);
v = zeros(1,time.stepNum);
U(:,1) = init.cellCenters;
v(1) = v_init;


for t = 2:nt
	t_step = t	
	Uprev = U(:,t-1);
	Vprev = v(t-1);

	for i = 1:10
		% Update Pennes Operator
		cutoff = CutoffFunction(Uprev,u_cutoff);
		scaling_cof = perfusion_scaling(Vprev,v_coma,v_star);
		pennes_op = DiagOperator(10^5 * dx,scaling_cof * pennes_cof.cellCenters,cutoff.vector);
		
		% Update Stiffness Operator
		stiff = StiffnessMatrixTime(dom,time,laplacian,boundary,pennes_op);

		% Update Source Vector
		F = tau * dom.cellWidths .* source.cellCenters + tau * boundary.source;
		F = [F; mu];
		F = F + [10^5 * dx .* U(:,t-1); v(t-1)];
		
		% Update ODE Vectors
		sense = SenseFunction(dom,Uprev,u_cutoff);
		pennes_vec = ColumnVector(dx,scaling_cof * pennes_cof.cellCenters,cutoff.vector);

		% Create Coupled Stiffness Matrix
		S = StiffnessMatrixCoupled(stiff,sense,pennes_vec,alpha).matrix;

		% Solve
		W = S \ F;
		Ucur = W(1:M);
		Vcur = W(M+1);

		% Update Iteration
		Uprev = Ucur;
		Vprev = Vcur;
	end

	U(:,t) = Ucur;
	v(t)   = Vcur;
end

U, v
