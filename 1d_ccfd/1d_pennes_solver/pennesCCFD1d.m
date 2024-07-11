function [U,xc,x,S] = pennesCCFD1d(domain,parameters,boundary,f)
%PENNESCCFD1D(domain,k,f,leftBC,leftBCtype,rightBC,rightBCtype) solves 1d 
%   stationary Pennes equation using cell-centered finite difference method
%
%   PENNESCCFD1D solves the equation -ku'' = f on the domain (0,L). 
%   Boundary conditions at x = 0 and x = L are given by leftBC and rightBC 
%   and may be Dirichlet type (indicated 'D') or Neumann type (indicated 
%   'N').
%
%   PENNESFE1D returns approximate solution U on domain x.
%
%   Author: Tyler Fara           Date: Sept 27, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   EXAMPLES
%   example 1
%       [U,xc] = pennesCCFD1d({1,10000},{1,1,1,0},{0,'D',0,'D'}); plot(xc,U);
%
%   example 2
%       f = @(x)(-2*exp(-x))
%       [U,x] = pennesCCFD1d(2,100,f,-2,'N',0,'D'); plot(x,U);
%
%   example 3 (remove extra line breaks from f when copy/pasting)
%      f = @(x)(-exp(1+(x.^2 - 2*x).^(-1)).*(x.^2 - 2*x).^(-4).*(2*x-2).^2 
%           - 2*exp(1+(x.^2-2*x).^(-1)).*(x.^2-2*x).^(-3).*(2*x-2).^2 
%           + 2*exp(1+(x.^2-2*x).^(-1)).*(x.^2-2*x).^(-2)); 
%      syms y(x); y(x) = piecewise(x<=0,0,0<x&x<2,f(x),2<=x,0);
%      [U,x] = pennesCCFD1d(1,100,y,0,'N',1,'D'); plot(x,U);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% VARIABLES
[k,c,theta_hat] = assignParameters(parameters);
perm  = Permeability(domain,k);
C	  = Coefficient(domain,c);
theta = Coefficient(domain,theta_hat);
tx    = Transmissibility(domain,perm);

% SETUP SU = Q
laplacian = DiscreteLaplacian(tx);

S = constructStiffnessMatrix(domain.cellWidths',tx.cellEdges',C.cellCenters');
Q = constructRHS(domain.cellWidths',domain.cellCenters,f,C.cellCenters',theta.cellCenters');

% BOUNDARY CONDITIONS
% error handling
checkBoundaryConditionFormat(boundary);

% update A and F with boundary conditions
S = updateStiffnessMatrix(tx.cellEdges',S,boundary);
Q = updateRHS(tx.cellEdges',Q,boundary);


% SOLVE
temp = S \ Q';
U = temp';

xc = domain.cellCenters';

end


%% DOMAIN FUNCTIONS

function [x,dx,xc] = constructDomain(domain)
%CONSTRUCTDOMAIN(domain) directs the construction of the domain

    x  = setGridPoints(domain);
    dx = setGridLengths(x);
    xc = setGridMidpoints(x,dx);

end


function x = setGridPoints(domain)
%SETGRIDPOINTS(domain) builds a vector of grid points by parsing input type

    % domain input is a cell array, i.e. {L,n} or {a,b,n}, where:
    %   L = right endpoint of domain, n = number of cells
    if isa(domain, "cell")
        if length(domain) == 2
            L = domain{1};
            n = domain{2};
            x = linspace(0,L,n+1);
        elseif length(domain) == 3
            a = domain{1};
            b = domain{2};
            n = domain{3};
            x = linspace(a,b,n+1);
        else
            error 'Incorrect domain input. Cell array should be format {L, n} or {a, b, n}.'
        end
        
    % domain input is a vector
    elseif isa(domain, 'double')
        if length(domain) == 1
            error 'Incorrect domain input. Vector input should have length at least two.'
        else
            x = domain;
        end
    end
end


function dx = setGridLengths(x)
%SETGRIDLENGTHS builds a vector of lengths of the cells of input x

    dx = diff(x);

end


function xc = setGridMidpoints(x,dx)
%SETGRIDMIDPOINTS(x) builds a vector of cell midpoints
    
    n  = length(x)-1;
    xc = x(1:n) + dx/2;

end


%% GENERAL FUNCTIONS

function tx = computeTransmissibility(dx,perm)
%COMPUTETRANSMISSIBILITY computes the transmissibility at each edge

    n    = length(dx);
    temp = zeros(n+1,1);
    temp(1:n) = perm./dx;
    tx   = zeros(n+1,1);

    for i = 2:n
        tx(i) = harmonicMean(temp(i-1),temp(i));
    end
    tx(1) = 2*temp(1);
    tx(n+1) = 2*tx(n);

end


function F = computeFunction(x,f)
%STORESOURCE evaluates f at x and stores the result in vector F.

    % parse class of f and evaluate f(x) to obtain F
    if isa(f,"double")
        F = f*ones(1,length(x));
    elseif isa(f,"function_handle")
		if nargin(f) == 1
        	F = f(x);
		elseif nargin(f) == 0
			F = f() * ones(1,length(x));
		end
    elseif isa(f,"symfun")
        F = f(x);
        F = double(F);
    else
        error 'Incorrect function type. Expected type is double, symfun, or function_handle.'
    end
end


function mean = harmonicMean(x,y)
%HARMONICMEAN(x,y) computes harmonic mean of x and y
    
    mean = 2/(1/x + 1/y);

end


%% SETTING UP PARAMETERS

function [k,c,theta,f] = assignParameters(parameters)

    k     = parameters{1};
    c     = parameters{2};
    theta = parameters{3};

end


%% STIFFNESS MATRIX FUNCTIONS

function S = constructStiffnessMatrix(dx,tx,C)
%CONSTRUCTSTIFFNESSMATRIX directs the construction of the stiffness 
%   matrix.

    n = length(tx)-1;
    A = sparse(n,n);
    
    for i = 2:n
        gl = i-1;
        gr = i;
    
        A(gl,gl) = A(gl,gl) + tx(i);
        A(gl,gr) = A(gl,gr) - tx(i); 
        A(gr,gl) = A(gr,gl) - tx(i);
        A(gr,gr) = A(gr,gr) + tx(i);
    end

    % B = sparse(diag(C).*dx);
    
    B = spdiags(C'.*dx,0,length(C),length(C));

    S = A + B;
end


%% RHS FUNCTIONS

function Q = constructRHS(dx,xc,f,C,theta)
%CONSTRUCTRHS directs the construction of the right hand side.

    F = computeFunction(xc,f);
    Q = F.*dx + theta.*C.*dx;
    
end


%% BOUNDARY CONDITION FUNCTIONS

function checkBoundaryConditionFormat(boundary)
%CHECKBOUNDARYCONDITIONFORMAT checks that the boundary condition types are 
%   formatted as 'D' or 'N'.
    
    % unpack boundary variables
	leftBC		= boundary{1};
	leftBCtype  = boundary{2};
	rightBC		= boundary{3};
	rightBCtype = boundary{4};

    % check boundary conditions
    if (leftBCtype == 'D' || leftBCtype == 'N' || leftBCtype == 'R') ...
        && (rightBCtype == 'D' || rightBCtype == 'N' || rightBCtype == 'R')
        ...
    else
        error 'Incorrect boundary condition types. Input ''D'' or ''N'' for boundary condition type.'
    end

	errmsg = 'Robin boundary condition should be of form q = alpha(u - uBar) and must be input as [alpha,uBar]';
	if leftBCtype == 'R' && length(leftBC) ~= 2 
		error errmsg
	elseif rightBCtype == 'R' && length(rightBC) ~=2
		error errmsg
	end

end


function checkWellPosedness(boundary)
%CHECKWELLPOSEDNESS checks that the problem is not ill-posed.

    % unpack boundary variables
    leftBCtype  = boundary{2};
    rightBCtype = boundary{4};

    % check well-posedness
    if (leftBCtype == 'N' && rightBCtype == 'N')
        error 'Ill-posed problem. Only one boundary may be type ''N''.'
    end
end


function S = updateStiffnessMatrix(tx,S,boundary)    
%UPDATESTIFFNESSMATRIX updates stiffness matrix S based on the 
%   boundary conditions.

    % unpack boundary variables
    leftBC      = boundary{1};
    leftBCtype  = boundary{2};
    rightBC     = boundary{3};
    rightBCtype = boundary{4};

    % store # of intervals
    n = length(tx)-1;

    % apply x = 0 boundary condition
    if leftBCtype == 'D'
        S(1,1) = S(1,1) + tx(1);
    elseif leftBCtype == 'N'
        ...
    elseif leftBCtype == 'R'
        S(1,1) = S(1,1) + tx(1)/(1 + tx(1)*(1/leftBC(1)));
    end
    
    % apply x = L boundary condition
    if rightBCtype == 'D'
        S(n,n)   = S(n,n) + tx(n+1);
    elseif rightBCtype == 'N'
        ...
    elseif rightBCtype == 'R'
        S(n,n) = S(n,n) + tx(n+1)/(1 + tx(n+1)*(1/rightBC(1)));
    end
end


function Q = updateRHS(tx,Q,boundary)
%UPDATERHS updates vector F, which represents the right-hand side, based 
%   on the boundary conditions.

    % unpack boundary variables
    leftBC      = boundary{1};
    leftBCtype  = boundary{2};
    rightBC     = boundary{3};
    rightBCtype = boundary{4};

    % store # of intervals
    n = length(tx)-1;

    % apply x = 0 boundary condition
    if leftBCtype == 'D'
        Q(1) = Q(1) + tx(1)*leftBC;
    elseif leftBCtype == 'N'
        Q(1) = Q(1) - leftBC;
    elseif leftBCtype == 'R'
        Q(1) = Q(1) + tx(1) / (1 + tx(1) * (1/leftBC(1))) * leftBC(2);
    end
    
    % apply x = L boundary condition
    if rightBCtype == 'D'
        Q(n) = Q(n) + tx(n+1)*rightBC;
    elseif rightBCtype == 'N'
        Q(n) = Q(n) - rightBC;
    elseif rightBCtype == 'R'
        Q(n) = Q(n) + tx(n+1) / (1 + tx(n+1)*(1/rightBC(1))) * rightBC(2);
    end
end
