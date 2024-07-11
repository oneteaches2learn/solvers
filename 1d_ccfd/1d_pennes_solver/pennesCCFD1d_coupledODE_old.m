function [U,uBody,xc,tGrid,S,handBoundaryIndArchive] = pennesCCFD1d_coupledODE3(domain,parameters,boundary,time,iter,nonlinearityflag)
%PENNESCCFD1DTIME(domain,parameters,boundary,time) computes a cell-centered
%   finite difference solution u on grid centers xc to the equation
%   u_t -k(x)u_xx + c(x)(u - u*) = f(x). 
%
%   Author: Tyler Fara           Date: Oct 27, 2022
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


% additional variables (to be handled differently later)
loop_verbosity = "on";
wait_verbosity = "on";
handBoundary = 0.05;
uAir      = -40;
uBodyInit = 37;
c_air     = 10^6;
c_hand    = 10^3;
c_B       = .5;
if nonlinearityflag == 1
	globalTempEffect = 1;
	localTempEffect  = 1;
	domainShrinking  = 1;
else
	globalTempEffect = 0;
	localTempEffect  = 0;
	domainShrinking  = 0;
end

gam = 10^2;

if localTempEffect == 1
	uCutoff = 20;
else 
	uCutoff = -40;
end

% VARIABLES
[L,M,x,dx,xc]      = constructDomain(domain);
[k,c,theta,f,por]  = assignParameters(parameters);
[dt,nt,tGrid,init] = constructTimeStepping(time);
[atol,maxIter]     = constructIteration(iter);
[leftCond,leftType,rightCond,rightType] = constructBoundary(boundary);

perm = computeFunction(k,xc);
Por  = computeFunction(por,xc);
Init = computeFunction(init,xc);
tx   = computeTransmissibility(dx,perm);

U     = zeros(nt+1,length(xc));
uBody = zeros(1,nt+1);
C     = zeros(1,length(xc));
uStar = zeros(1,length(xc));



% SOLVER
% Apply initial condition
U(1,:)    = Init;
handBoundaryInd = get_handBoundary_index(xc,handBoundary,U(1,:),uCutoff);
uBody(1)  = uBodyInit;
uArterial = uBody(1);
uStar     = computeUStar(uStar, xc, handBoundary, uAir, uBody(1));
C         = computeNonlinearPennesCoefficient(...
				C,xc,U(1,:),uCutoff,uBody(1),handBoundary,c_air,c_hand); 

% Ode setup
mu = 2; lambda = 0;
tRange_ODE = [0,nt*dt];
y0_ODE     = uBody(1);
uSense     = get_uSense(handBoundaryInd,U(1,:));
if globalTempEffect == 1
	ode = @(t,y) my_nonlinear_ODE(t,y,mu,lambda,c_B,uSense);
else
	ode = @(t,y) my_linear_ODE(t,y,mu,lambda,c_B,uSense);
end
sol_ODE    = ode45(ode,tRange_ODE,y0_ODE);
uBody(1)   = deval(sol_ODE,0);

% Build base stiffness matrix 
S_base = constructStiffnessMatrix(dx,tx,Por,dt);

handBoundaryIndArchive = zeros(nt,1);
handBoundaryIndArchive(1) = handBoundaryInd;

% loop over time
if wait_verbosity == 'on'; wait = waitbar(0,'Please wait...'); end
if loop_verbosity == 'on'; fprintf('\n'); end
for i = 2:nt+1

    % Set up iteration
    t = (i-1) * dt;
	tRange_ODE = [t-dt,t+dt];
    iter = 1;
    Uiter(1,:) = U(i-1,:);
	y0_ODE = uBody(i-1);
    if loop_verbosity == "on"
        fprintf('     Computing timestep %d of %d...',i-1,nt);
    end

	if domainShrinking == 1
		handBoundaryInd = get_handBoundary_index(xc,handBoundary,Uiter(1,:),uCutoff);
	end

	handBoundaryIndArchive(i) = handBoundaryInd;

    % Iterate to resolve nonlinearity
    while 1

        % UPDATE NONLINEARITY
		% Update uBody
		uSense    = get_uSense(handBoundaryInd,Uiter(1,:));
		if globalTempEffect == 1
			ode = @(t,y) my_nonlinear_ODE(t,y,mu,lambda,c_B,uSense);
		else
			ode = @(t,y) my_linear_ODE(t,y,mu,lambda,c_B,uSense);
		end
		sol_ODE   = ode45(ode,tRange_ODE,y0_ODE);
		uBodyIter = deval(sol_ODE,t);
		if uBodyIter > 37
			uBodyIter = 37;
		end
		if uBodyIter < 32 && globalTempEffect == 1;
			uSense = 0;
			ode = @(t,y) my_nonlinear_ODE(t,y,mu,lambda,c_B,uSense);
			sol_ODE   = ode45(ode,tRange_ODE,y0_ODE);
			uBodyIter = deval(sol_ODE,t);
		end
		uArterial = uBodyIter;
		
		% Update variables
        F = computeFunction(f,xc,t,Uiter(1,:));
 		uStar = computeUStar(uStar, xc, handBoundary, uAir, uArterial);
		boundary = {leftCond,leftType,[gam,uBodyIter],rightType};
		if globalTempEffect == 1 || localTempEffect == 1
  			C = computeNonlinearPennesCoefficient(...
					C,xc,Uiter(1,:),uCutoff,uBodyIter,handBoundary,c_air,c_hand); 
		end

		% Update stiffness matrix
		S = updateStiffnessMatrixNonlinearity(S_base,dt,dx,C);
		S = updateStiffnessMatrixBoundary(tx,S,boundary,dt);

        % Update RHS
		G = constructRHS(dx,xc,t,dt,F,C,uStar,Por,U(i-1,:));
        G = updateRHSBoundary(tx,G,boundary,dt);


        % SOLVE
        temp = S \ G';
        Uiter(2,:) = temp';


        % CHECK STOPPING CRITERION
		[U,Uiter,uBody,iterCase,iter] = checkStoppingCriterion( ...
						U,Uiter,uBody,uBodyIter,L,M,atol,maxIter,loop_verbosity,i,iter);

		if iterCase == 0;     break;	% convergence acheived
		elseif iterCase == 1; break;	% no convergence before maxIter
		elseif iterCase == 2; ... 		% continue iterating
		end

	end
    if wait_verbosity == 'on'; waitbar(i/nt, wait, 'Please wait...'); end
end
if wait_verbosity == 'on'; close(wait); end

end 
% End of main program



%% ODE FUNCTIONS

function dydt = my_nonlinear_ODE(t,y,mu,lambda,c_B,uSense)

	if uSense < 20
		M = 0;
	elseif y > uSense
		M = y - uSense;
	elseif y <= uSense
		M = 0;
	end 
	dydt = (mu - lambda) - c_B * M;
%	dydt = (mu - lambda) - c_B * (y - uSense);

end

function dydt = my_linear_ODE(t,y,mu,lambda,c_B,uSense)

	dydt = (mu - lambda) - c_B * (y - uSense);

end

%% DOMAIN FUNCTIONS

function [L,M,x,dx,xc] = constructDomain(domain)
%CONSTRUCTDOMAIN(domain) directs the construction of the domain

    L  = domain{1};
    M  = domain{2};
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


%% SETTING UP PARAMETERS

function [k,c,uBody,f,por] = assignParameters(parameters)

    k     = parameters{1};
    c     = parameters{2};
    uBody = parameters{3};
    f     = parameters{4};
    por   = parameters{5};

end

function [atol,maxIter] = constructIteration(iter)

    atol    = iter{1};
    maxIter = iter{2};

end

function C = computeNonlinearPennesCoefficient(C,xc,U,uCutoff,uBody,handBoundary,c_air,c_hand)
%UPDATENONLINEARPENNESCOEFFICIENT(C,XC,U,UCUTOFF,HANDBOUNDARY,C_AIR,C_HAND)
%Rcomputes the Pennes coefficient C by comparing XC to HANDBOUNDARY and assigning
%	C = C_AIR for XC outside the hand and C = C_HAND or C = 0 for XC inside the hand.
%	For XC inside the hand, the determination of C in each element depends on whether
%	U > UCUTOFF or U < UCUTOFF.
 
	for i = 1:length(xc)
		if xc(i) < handBoundary 
			C(i) = c_air;
		elseif xc(i) >= handBoundary && uBody < 32
			C(i) = 0;
		elseif xc(i) >= handBoundary && U(i) >= uCutoff
			C(i) = c_hand * (37 - uBody) / 5;
		else
			C(i) = 0;
		end
	end

end

function uStar = computeUStar(uStar, xc, handBoundary, uAir, uArterial);
%COMPUTEUSTAR(USTAR,XC,HANDBOUNDARY,UAIR,UBODY) computes USTAR by comparing XC to
%	HANDBOUNDARY and assigning USTAR = UAIR for XC outside the hand and
%	USTAR = UARTERIAL for XC inside the hand.

 	for i = 1:length(xc)
	    if xc(i) < handBoundary
        	uStar(i) = uAir;
        end
        if xc(i) >= handBoundary
    	    uStar(i) = uArterial;
        end
	end

end


function uSense = get_uSense(handBoundaryInd, U)

	if handBoundaryInd < length(U)
		uSense = mean(U(handBoundaryInd+1:length(U)));
	else
		uSense = 0;
	end

end


function ind = get_handBoundary_index(xc,handBoundary,U,uCutoff)
%GET_HANDBOUNDARY_INDEX(XC,HANDBOUNDARY) returns the index of the first cell
%	of the hand. 
	
	% Find last index before the handBoundary
	for i = 1:length(xc)
		if xc(i) < handBoundary
			ind = i;
		elseif xc(i) >= handBoundary && U(i) < uCutoff
			ind = i;
		elseif xc(i) >= handBoundary && U(i) >= uCutoff			
			break;
		end
	end
end



%% TIME STEPPING

function [dt,nt,tGrid,init] = constructTimeStepping(time)

    dt    = time{1};
    nt    = time{2};
    init  = time{3};
    tGrid = linspace(0,nt*dt,nt+1);

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


function F = computeFunction(f,varargin)
%COMPUTEFUNCTION evaluates f(x) and stores result as vector F or evalutes
%   f(x,t) and stores result as matrix F.
    if nargin == 2
        x = varargin{1};
        F = evaluateOneInput(f,x);
    elseif nargin == 3
        x = varargin{1};
        t = varargin{2};
        F = evaluateTwoInputs(f,x,t);
    elseif nargin == 4;
        x = varargin{1};
        t = varargin{2};
        U = varargin{3};
        F = evaluateThreeInputs(f,x,t,U);
    else
        error 'Too many arguments passed. Expected two to four arguments.'
    end
end


function F = evaluateOneInput(f,x)
%EVALUATEONEINPUT evaluates f at x and stores the result in vector F.
    % parse class of f and evaluate f(x) to obtain F

    if isa(f,"double") && length(f) == 1
        F = f*ones(1,length(x));
    elseif isa(f,"function_handle")
        F = f(x);
    elseif isa(f,"symfun")
        F = f(x);
        F = double(F);
	elseif (isvector(f) && size(f) == size(x))
		F = f;
    else
        error 'Incorrect function type. Expected type is double, symfun, or function_handle.'
    end
end


function F = evaluateTwoInputs(f,x,t)
%EVALUATETWOINPUTS parses function class, evaluates f(x,t), and stores 
%   result in matrix F.
    % parse class of f and evaluate f(x) to obtain F

	if isa(f,"double") && length(f) == 1
        F = f*ones(1,length(x));
    elseif isa(f,"function_handle")
        F = f(x,t);
    elseif isa(f,"symfun")
        F = f(x,t);
        F = double(F);
	elseif isvector(f)
		F = f;
    else
        error 'Incorrect function type. Expected type is double, symfun, or function_handle.'
    end
end


function F = evaluateThreeInputs(f,x,t,U)
%EVALUATETWOINPUTS parses function class, evaluates f(x,t), and stores 
%   result in matrix F.

    % parse class of f and evaluate f(x) to obtain F
	if isa(f,"double") && length(f) == 1
        F = f*ones(1,length(x));
    elseif isa(f,"function_handle")
        F = f(x,t,U);
    elseif isa(f,"symfun")
        F = f(x,t,U);
        F = double(F);
    elseif isvector(f);
		F = f;
	else
        error 'Incorrect function type. Expected type is double, symfun, or function_handle.'
    end
end



function mean = harmonicMean(x,y)
%HARMONICMEAN(x,y) computes harmonic mean of x and y
    
    mean = 2/(1/x + 1/y);

end

function [U,Uiter,uBody,iterCase,iter] = checkStoppingCriterion(...
								U,Uiter,uBody,uBodyIter,L,M,atol,maxIter,loop_verbosity,i,iter)

    residual = sqrt(sum((L / (M + 1)) .* (Uiter(2,:) - Uiter(1,:)).^2));

    % If convergence was acheived
    if residual < atol
		iterCase = 0;
		U(i,:)   = Uiter(2,:);
		uBody(i) = uBodyIter;
		if loop_verbosity == "on"
        	fprintf(' converged in %g iterations.\n',iter);
   		end
        
    % Else if too many iterations have occurred
    elseif iter >= maxIter 
		iterCase = 1;
		uBody(i) = uBodyIter;
        fprintf('   No convergence, timestep %g. \n\n',i-1); 
        
    % Else prepare for next iteration    
    else
		iterCase = 2;
        iter = iter+1; 
        Uiter(1,:) = Uiter(2,:); 
    end

end


%% STIFFNESS MATRIX FUNCTIONS

function S = constructStiffnessMatrix(dx,tx,Por,dt)
%CONSTRUCTSTIFFNESSMATRIX directs the construction of the stiffness 
%   matrix.

    S1 = buildBaseMatrix(tx);
    S3 = sparse(length(Por),length(Por));

    for i = 1:length(Por)
        S3(i,i) = Por(i) * dx(i);
    end

    S = dt*S1 + S3;

end


function A = buildBaseMatrix(tx)
%BUILDBASEMATRIX builds the tridiagonal matrix using transmissibility. 

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
end


function S_updated = updateStiffnessMatrixNonlinearity(S,dt,dx,C);
%UPDATENONLINEARSTIFFNESSMATRIX(S,DT,DX,C) adds the Pennes coefficient tensor
%	to the stiffness matrix tensor.
		
		S_updated = S;  
        for i = 1:length(C)
            S_updated(i,i) = S_updated(i,i) + dt * C(i) * dx(i);
        end

end



%% RHS FUNCTIONS

function Q = constructRHS(dx,xc,t,dt,f,C,theta,Por,nsol)
%CONSTRUCTRHS directs the construction of the right hand side.

    F = computeFunction(f,xc,t);
    Q = dt*F.*dx + dt*theta.*C.*dx + Por.*dx.*nsol;
    
end


%% BOUNDARY CONDITION FUNCTIONS

function S = updateStiffnessMatrixBoundary(tx,S,boundary,dt)    
%UPDATESTIFFNESSMATRIX updates stiffness matrix S based on the 
%   boundary conditions.


    % error handling for boundary conditions
    checkBoundaryConditionFormat(boundary);
    checkWellPosedness(boundary);

    % unpack boundary variables
    leftBC      = boundary{1};
    leftBCtype  = boundary{2};
    rightBC     = boundary{3};
    rightBCtype = boundary{4};

    % store # of intervals
    n = length(tx)-1;

    % apply x = 0 boundary condition
    if leftBCtype == 'D'
        S(1,1) = S(1,1) + dt*tx(1);
    elseif leftBCtype == 'N'
        ...
    elseif leftBCtype == 'R'
        S(1,1) = S(1,1) + dt*leftBC(1);
    end
    
    % apply x = L boundary condition
    if rightBCtype == 'D'
        S(n,n)   = S(n,n) + dt*tx(n+1);
    elseif rightBCtype == 'N'
        ...
    elseif rightBCtype == 'R'
        S(n,n) = S(n,n) + dt*rightBC(1);
    end
end


function [leftCond,leftType,rightCond,rightType] = constructBoundary(boundary)

	leftCond  = boundary{1};
	leftType  = boundary{2};
	rightCond = boundary{3};
	rightType = boundary{4};

end


function checkBoundaryConditionFormat(boundary)
%CHECKBOUNDARYCONDITIONFORMAT checks that the boundary condition types are 
%   formatted as 'D' or 'N'.
    
    % unpack boundary variables
    leftBC      = boundary{1};
    leftBCtype  = boundary{2};
    rightBC     = boundary{3};
    rightBCtype = boundary{4};

    % check boundary conditions
    if (leftBCtype == 'D' || leftBCtype == 'N') ...
        && (rightBCtype == 'D' || rightBCtype == 'N')
        ...
    else
    %    error 'Incorrect boundary condition types. Input ''D'' or ''N'' for boundary condition type.'
    end
end


function checkWellPosedness(boundary)
%CHECKWELLPOSEDNESS checks that the problem is not ill-posed.

    % unpack boundary variables
    leftBCtype  = boundary{2};
    rightBCtype = boundary{4};

    %{
    % check well-posedness
    if (leftBCtype == 'N' && rightBCtype == 'N')
        error 'Ill-posed problem. Only one boundary may be type ''N''.'
    end
    %}
end


function Q = updateRHSBoundary(tx,Q,boundary,dt)
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
        Q(1) = Q(1) + dt*tx(1)*leftBC;
    elseif leftBCtype == 'N'
        Q(1) = Q(1) - dt*leftBC;
    elseif leftBCtype == 'R'
        Q(1) = Q(1) + dt*leftBC(1)*rightBC(2);
    end
    
    % apply x = L boundary condition
    if rightBCtype == 'D'
        Q(n) = Q(n) + dt*tx(n+1)*rightBC;
    elseif rightBCtype == 'N'
        Q(n) = Q(n) - dt*rightBC;
    elseif rightBCtype == 'R'
        Q(n) = Q(n) + dt*rightBC(1)*rightBC(2);
    end
end
