function [U,xc,tGrid] = pennesCCFD1dtime(domain,parameters,boundary,time)
%PENNESCCFD1DTIME(domain,parameters,boundary,time) computes a cell-centered
%   finite difference solution u on grid centers xc to the equation
%   u_t - k(x)u_xx + c(x)(u - u*) = f(x). 
%
%   Author: Tyler Fara           Date: Feb 29, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUTS
%       domain:     {L,M}
%       parameters: {k,c,theta,f,por} 
%       boundary:   {leftCond,leftType,rightCond,rightType}
%       time:       {dt,nt,init} 
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


% VARIABLES
% unpack inputs
[x,dx,xc]          = constructDomain(domain);
[k,c,theta,f,por]  = assignParameters(parameters);
[dt,nt,tGrid,init] = constructTimeStepping(time);

% compute variables
perm      = computeFunction(k,xc);
C         = computeFunction(c,xc);
Por       = computeFunction(por,xc);
Init      = computeFunction(init,xc);
theta     = computeFunction(theta,xc);
tx        = computeTransmissibility(dx,perm);


% BUILD STIFFNESS MATRIX
S = constructStiffnessMatrix(dx,tx,C,Por,dt);
S = updateStiffnessMatrix(tx,S,boundary,dt);


% SOLVER
% apply initial condition
U = zeros(nt+1,length(xc));
U(1,:) = Init;

% loop over time
%wait = waitbar(0,'Please wait...');
for i = 2:nt+1

    t = i*dt;

	% create RHS vector
    Q = constructRHS(dx,xc,t,dt,f,C,theta,Por,U(i-1,:));
    Q = updateRHS(tx,Q,boundary,dt,t);

    % solve
    temp = S \ Q';
    U(i,:) = temp';
    
%    waitbar(i/nt, wait, 'Please wait...');
end

%close(wait);
end




















