% DEMO1_DRIVER 
%	On rectangle Omega = xLim x yLim, let u = uTrue satisfy  
%
%		 -Delta u + u = f, in Omega
%		            u = u_D on Omega_D
%		-grad u dot n = u_N on Omega_N
%		-grad u dot n = alpha * (u - u_R) on Omega_R
%
%	where uTrue is specified by the user. demo1_driver manufactures the
%	necessary boundary conditions and source terms in order to generate uTrue
%	from the GalerkinPoisson2d_solver. 
%
%	author: Tyler Fara						date: June 27, 2024
%-----------------------------------------------------------------------------%
% NOTES
%	(1) Omega is triangulated using triangles with diameter less than h =
%	base^-p.
%
%	(2) If variable mms_test = 1, then demo1_driver runs an mms test with p =
%	{1, 2, 3, 4}. In this case, the user's choice of p is disregarded. If
%	mms_test = 0, then demo1_driver manufactures a single solution for h =
%	base^p. 
%
%	(3) Acceptable choices of bTypes are 'D' for Dirichlet, 'N' for Neumann, or
%	'R' for Robin. The choices corespond to the southern edge, the eastern
%	edge, the northern edge, and the western edge respectively. All
%	combinations of boundary conditions should work.
%
%	(4) If running mms test, you should observe order 2 convergence. The result
%	will be stored as an object named mms, and the data stored in mms will be
%	printed to the command line. Look for the numbers next to the property
%	'orders' to converge toward the number 2.
%
%	(5) You should be able to specify any continuous uTrue. uTrue is specified
%	as a symbolic function with variables x(1), representing the x direction,
%	and x(2), representing the y direction.
%	
%	(6) Try setting mms_test = 0. Then try prob.plot. Or try prob.domain.plot.
%	Or prob.domain.plot(NodeLabels="on",ElementLabels="on")
%
%-----------------------------------------------------------------------------%


clear all; x = sym('x',[1 2],'real');
% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain bounds
xLim = [0 1];
yLim = [0 1];

% generate mesh
p = 2;
base = 2;

% specify BCs
bTypes = {'R' 'R' 'R' 'R'};

% specify desured result
uTrue = sin(pi/2*x(1))*sin(pi/2*x(2));

mms_test = 1;


% MMS TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ensure functions have correct variables
uTrue = symfun(uTrue,x);

% assemble inputs
dom   = Domain2d(xLim,yLim,p,base);
bound = Boundary2d(bTypes,{});
cofs  = [];  % <~~coefficients not yet implemented

% run mms test
if mms_test == 1
	mms = GalerkinPoisson2d_mms(dom,bound,cofs,uTrue,errType="L2")

% run demo test
elseif mms_test == 0
	mms = GalerkinPoisson2d_mms(dom,bound,cofs,uTrue,demo=p);
	prob = mms.problems{1};
end

