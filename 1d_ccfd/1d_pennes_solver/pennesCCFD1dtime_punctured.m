function [U,xc,tGrid] = pennesCCFD1dtime_punctured(domain,parameters,boundary,time,punctures)
%PENNESCCFD1DTIME_PUNCTURED(domain,parameters,boundary,time) computes a cell-centered
%   finite difference solution u on grid centers xc to the equation
%   u_t - k(x)u_xx + c(x)(u - u*) = f(x). This version allows for punctures in the 
%	domain, i.e. of the form [0,a] U [b,L]
%
%   Author: Tyler Fara           Date: Feb 29, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUTS
%       domain:     {L,M,punc1}
%       parameters: {k,c,theta,f,por} 
%       boundary:   {leftCond,leftType,rightCond,rightType,punc1a,punc1aType,punc1b,punc1bType}
%       time:       {dt,nt,init} 
%		varargin:	one or more Puncture.m objects
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   EXAMPLES
%   example 1
%       L = 1; M = 100; dt = 1; nt = 10; 
%       [U,xc,tGrid] = pennesCCFD1dtime({L,M},{1,1,1,0,1},{0,'D',0,'D'},{dt,nt,0}); 
%       surfplotterCCFD(U,L,dt,10);
%
%   example 2
%       init = @(x)sin(pi*x); L = 1; dt = .1;
%       [U,xc] = pennesCCFD1dtime({L,100},{1,1,1,0,1},{0,'D',0,'D'},{dt,10,init});
%       surfplotterCCFD(U,L,dt,10);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Unpack inputs
a = domain{1};
b = domain{2};
N = domain{3};
leftBC    = boundary{1};
leftType  = boundary{2};
rightBC   = boundary{3};
rightType = boundary{4};
dt = time{1};
nt = time{2};

% compute necessary variables
L = (b-a);
h = L/N;
subdom_num = length(punctures) + 1;

% set up subdomain endpoints and BCs
[BC_L,BCtype_L,BC_R,BCtype_R] = set_subdomain_BCs(boundary,punctures);
[endpoints_L,endpoints_R]     = set_subdomain_endpoints(a,b,punctures);

% Compute solution on each subdomain
U = []; xc = [];
for i = 1:subdom_num

	% Divide subinterval into cells of length at most h
	x_L = endpoints_L(i);
	x_R = endpoints_R(i);
	N_i = ceil((x_R - x_L) / h); % number of cells, ensuring that cell width < h
	sub = linspace(x_L,x_R,N_i+1);

	% Solve on subdomain
	domain = {endpoints_L(i),endpoints_R(i),N_i};
	boundary = {BC_L{i},BCtype_L{i},BC_R{i},BCtype_R{i}};
	[U_sub,xc_sub,tGrid] = pennesCCFD1dtime(domain,parameters,boundary,time);
	
	% Assemble result
	if i < subdom_num, xc = [xc xc_sub NaN];
	else xc = [xc xc_sub]; end
	
	if i < subdom_num, U = [U U_sub NaN(nt+1,1)];
	else U = [U U_sub]; end

end

end
