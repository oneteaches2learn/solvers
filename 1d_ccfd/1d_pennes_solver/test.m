clear variables

% DATA
% Coefficients and Source
perm_fun   = @(x)(1 + 0*x);
pennes_fun = @(x)(1 + 0*x);
source_fun = @(x)(1 + 0*x);

% Time Data
tau = 0.1;
nt  = 10;
init_fun = @(x)(sin(pi*x));

% Domain Data
a = 0; b = 1; M = 11;
cellEdges = linspace(a,b,M);

% Boundary Data
u_a = 0;
u_b = 0;

% ODE data
alpha = 1;
mu = 0;
v_init = 10;


% SOLVER
% Construct domain and time-stepping
dom  = Domain(cellEdges);
time = TimeStepping(tau,nt); 

% Construct coefficients
init   = InitialCondition(dom,init_fun);
perm   = Permeability(dom,perm_fun);
source = Source(dom,perm_fun);
pennes = Coefficient(dom,pennes_fun);
trans  = Transmissibility(dom,perm);
laplacian = DiscreteLaplacian(trans);

% Construct Boundary Conditions
left     = DirichletBoundary(u_a,"L",trans);
right    = DirichletBoundary(u_b,"R",trans);
boundary = BoundaryConditions(left,right,dom);

M  = dom.cellNum;
dx = dom.cellWidths;
tau = time.stepSize;


% PENNES SOLVER
% Initialize data storage
U = zeros(dom.cellNum,time.stepNum);
v = zeros(1,time.stepNum);
U(:,1) = init.cellCenters;
v(1) = v_init;

% Create Stiffness Matrix and Source Vector
stiff = StiffnessMatrixTime(dom,time,laplacian,boundary,DiagOperator(dx,pennes.cellCenters));
S = stiff.matrix;

F = tau * dom.cellWidths .* source.cellCenters + tau * boundary.source;

for i = 1:time.stepNum
	nt = i+1;
	G = F + dx .* U(:,nt-1) + 10 * pennes.cellCenters .* dx;
	U(:,nt) = S \ G;
end

U


% COUPLED ODE SOLVER
% Initialize data storage
U = zeros(dom.cellNum,time.stepNum);
v = zeros(1,time.stepNum);
U(:,1) = init.cellCenters;
v(1) = v_init;

% Extend Stiffness Matrix
S_bar = zeros(M+1,M+1);
S_bar(1:M,1:M) = S;
S_bar(M+1,1:M) = -(alpha / dom.length) * dx';
S_bar(1:M,M+1) = -dx .* pennes.cellCenters;
S_bar(M+1,M+1) = alpha + 1;

% Extend Source Vector
F_bar = zeros(M+1,1);
F_bar(1:M) = F;
F_bar(M+1) = mu;

for i = 1:time.stepNum

	nt = i+1;
	G_bar = F_bar + [dx .* U(:,nt-1); v(nt-1)];
	temp = S_bar \ G_bar;

	U(:,nt) = temp(1:M);
	v(nt)   = temp(M+1);

end

U, v
