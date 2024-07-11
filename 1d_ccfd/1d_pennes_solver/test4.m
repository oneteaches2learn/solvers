clear variables

% DATA
% Coefficients and Source
perm_fun   = @(x)(1 + 0*x);
pennes_fun = @(x)(1 + 0*x);
source_fun = @(x)(0 + 0*x);

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


% SOLVER
% Construct domain and time-stepping
dom  = Domain(cellEdges);
time = TimeStepping(tau,nt); 

% Construct coefficients
init   = InitialCondition(dom,init_fun);
perm   = Permeability(dom,perm_fun);
source = Source(dom,perm_fun);
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
U(:,1) = init.cellCenters;

% Create Stiffness Matrix and Source Vector
stiff = StiffnessMatrixTime(dom,time,laplacian,boundary);
S = stiff.matrix;

F = tau * dom.cellWidths .* source.cellCenters + tau * boundary.source;

for i = 1:time.stepNum
	nt = i+1;
	G = F + dx .* U(:,nt-1);
	U(:,nt) = S \ G;
end

U
