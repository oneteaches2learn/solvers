function [U,xc,tGrid] = pennesCCFD1dTimeNonlinear(domain,parameters,boundary,time,iterParm)
%PENNESCCFD1DTIMENONLINEAR(domain,parameters,boundary,time) computes a
%	cell-centered finite difference solution u on grid centers xc to the equation
%   
%			por*u_t - k(x)u_xx + c(x)(u - uStar) = f(x). 
%
%   Author: Tyler Fara           Date: April 2, 2024
%------------------------------------------------------------------------------%
% Inputs
%	domain:     {L,M}
%	parameters: {k,c,theta,f,por} 
% 	boundary:   {leftCond,leftType,rightCond,rightType}
%	time:       {dt,nt,init} 
%	iterParm:	{atol,maxIter,verbosity}
%
% Outputs
%	U
%	xc
%	tGrid
%------------------------------------------------------------------------------%
% Notes
%	(1) As of April 3, 2024, this solver is semilinear, i.e. it can handle
%	nonlinearities in any coefficient except permeability, k. I plan to extend
%	it to a quasilinear solver. However, for now, semilinear is all I need. To
%	the user, k = k(x,t) will work; but k = k(x,t,u) may not work. 
%
%	(2) As of April 3, 2024, whether to use column vectors or row vectors is a
%	problem. I'd like to switch to the use of column vectors for all vectors
%	because this makes more intuitive sense, especially when it comes to
%	matrix-vector products. However, pennesCCFD1dtime, pennesCCFD1d, and
%	pennesCCFD1dtime_coupledODE all use row vectors, and several functions that
%	generate these row vectors are shared. So if I update those functions to
%	create column vectors so that this present nonlinear solver will work, that
%	will necessarily mean breaking the other solvers. Thus, all solvers need to
%	be updated together. That is something I'd like to do, but it is a problem
%	for another day!
%------------------------------------------------------------------------------%


	% unpack inputs
	[x,dx,xc]           = constructDomain(domain);
	[k,c,uStar,f,por]   = assignParameters(parameters);
	[dt,nt,tGrid,init]  = constructTimeStepping(time);

	% make xc and dx into column vectors...for now...
	dx = dx';
	xc = xc';

	% initialize storage
	U = zeros(nt+1,length(xc));   % stored as rows...for now...
	Uiter = zeros(2,length(xc));  % stored as rows...for now...

	% compute initial condition
	Init   = computeFunction(init,xc); % <--this is a column...
	U(1,:) = Init;					   % <--but is stored here as a row...

	% loop over time
	for i = 2:nt+1

		% set up iteration
		t = i*dt;
		Uiter(1,:) = U(i-1,:);
		iter = 0;

		% Iterate to resolve nonlinearity
		while 1

			% compute coefficients
			K     = computeFunction(k,xc,t,Uiter(1,:)');
			C     = computeFunction(c,xc,t,Uiter(1,:)');
			Por   = computeFunction(por,xc,t,Uiter(1,:)');
			UStar = computeFunction(uStar,xc,t,Uiter(1,:)');
			F     = computeFunction(f,xc,t,Uiter(1,:)');
			tx    = computeTransmissibility(dx,K);
			Q_bc  = computeRHSbc(tx,boundary,t);

			% compute tensors 
			M_k  = computeStiffnessMatrix(tx);
			A_h  = computeMassMatrix(Por,dx);
			M_c  = computeMassMatrix(C,dx);
			M_bc = computeBCMatrix(tx,boundary);

			% construct LHS
			S = A_h + dt*M_k + dt*M_c + dt*M_bc;

			% construct RHS
			Q = A_h*U(i-1,:)' + dt*M_c*UStar + dt*F.*dx + dt*Q_bc;
			dt*M_c*UStar

			% SOLVE
			temp = S \ Q;
			Uiter(2,:) = temp';

			% CHECK STOPPING CRITERION
			L = x(length(x)) - x(1);
			M = x(2) - x(1);
			[U,Uiter,iterCase,iter] = checkStoppingCriterionNonlinear(...
					U,Uiter,dx,iterParm,i,iter);

			if iterCase == 0;     break;	% convergence acheived
			elseif iterCase == 1; break;	% no convergence before maxIter
			elseif iterCase == 2; ... 		% continue iterating
			end
		end
	end
end

