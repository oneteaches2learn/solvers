clear all;
x = sym('x',[1 2],'real'); syms t;
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain bounds
xLim = [0 1];
yLim = [0 1];

% mesh parameters
p = 7;
base = 2;

% time parameters
N_t = 1001;

% specify function
f = sin(pi/2 * x(1)) * sin(pi/2 * x(2));


% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct domain

fprintf('Running Quadrature Test\n')
fprintf(' Contructing Domain:'), tic
	dom = GalerkinHeat2d_assembler.assembleDomainGeometry(xLim,yLim);
	dom = GalerkinHeat2d_assembler.assembleMesh(dom,p,base); 
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% check function
fprintf(' Assembling Parameters:'), tic
	f = symfun(f,x);
	f = matlabFunction(f);
executionTime = toc; 
fprintf(' %f s\n',executionTime)

% compute quadrature
fprintf(' Computing Quadrature:')
tic;
quad = parallel_quad_ver2(dom,f,N_t);
executionTime = toc; 
fprintf(' %f s\n',executionTime)
fprintf(' nElem: %d \n',dom.nElem)



