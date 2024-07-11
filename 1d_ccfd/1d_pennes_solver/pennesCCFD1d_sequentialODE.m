function [U,xc,tGrid,V] = pennesCCFD1d_sequentialODE(domain,parameters,boundary,time,iterParm,odeParm)
%PENNESCCFD1D_SEQUENTIALODE(domain,parameters,boundary,time) computes a
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
	[k,gam,uStar,f,por] = assignParameters(parameters);
	[dt,nt,tGrid,init]  = constructTimeStepping(time);
	[v_0,fB,gamB,S,vLower,vUpper] = assignODEparm(odeParm);
	
	% make xc and dx into column vectors...for now...
	dx = dx';
	xc = xc';

	% initialize storage
	U = zeros(nt+1,length(xc));   % stored as rows...for now...
	V = zeros(nt+1,1);
	Uiter = zeros(2,length(xc));  % stored as rows...for now...
	Viter = zeros(2,1);

	% compute initial condition
	Init   = computeFunction(init,xc); % <--this is a column...
	U(1,:) = Init;					   % <--but is stored here as a row...
	V(1) = v_0;

	% loop over time
	for i = 2:nt+1

		% set up iteration
		t = i*dt;
		iter = 1;
		Uiter(1,:) = U(i-1,:);
		Viter(1) = V(i-1);

		verbosity = 0;
		if verbosity == 1;
		if max(dx) > 0.05
			fprintf('  TimeStep %d\n',i)
		end
		end

		% Iterate to resolve nonlinearity
		while 1

			% compute body temperature
			Viter(2) = (dt * fB(i) + dt * gamB(i) * S(Uiter(1,:),dx) + V(i-1)) / (1 + dt * gamB(i));	 % Discrete uAvg
			%Viter(2) = (dt * fB(t) + dt * gamB * S(Uiter(1,:),dx) + V(i-1)) / (1 + dt * gamB);   % Integral uAvg

			% compute coefficients
			K     = computeFunction(k,xc,t,Uiter(1,:)',Viter(2));
			Gam   = computeFunction(gam,xc,t,Uiter(1,:)',Viter(2));
			Por   = computeFunction(por,xc,t,Uiter(1,:)',Viter(2));
			UStar = computeFunction(uStar,xc,t,Uiter(1,:)',Viter(2));
			F     = computeFunction(f,xc,t,Uiter(1,:)',Viter(2));
			tx    = computeTransmissibility(dx,K);
			Q_bc  = computeRHSbc(tx,boundary,t,Uiter(1,:),Viter(2));

			% compute tensors 
			A_h   = computeStiffnessMatrix(tx);
			M_c   = computeMassMatrix(Por,dx);
			M_gam = computeMassMatrix(Gam,dx);
			M_bc  = computeBCMatrix(tx,boundary,t,Uiter(1,:)',Viter(2));  

			% construct LHS
			A = M_c + dt*A_h + dt*M_gam + dt*M_bc;

			% construct RHS
			Q = M_c*U(i-1,:)' + dt*M_gam*UStar + dt*F.*dx + dt*Q_bc;

			% SOLVE
			temp = A \ Q;
			Uiter(2,:) = temp';

			% Check variables
			if verbosity == 1
			if max(dx) > 0.05
				fprintf('    Iteration %d\n',iter);
				fprintf('      Viter(2):                 %d\n',Viter(2));
				fprintf('      S(Uiter(1,:)):            %d\n',S(Uiter(1,:),dx)); 
				fprintf('      dt*gamB*S(Uiter(1,:),dx): %d\n',dt*gamB*S(Uiter(1,:),dx)); 
				fprintf('      dt*fB(t):                 %d\n',dt*fB(t)); 
			end
			end

			% CHECK STOPPING CRITERION
			L = x(length(x)) - x(1);
			M = x(2) - x(1);
			[U,Uiter,iterCase,iter] = checkStoppingCriterionNonlinear(...
					U,Uiter,dx,iterParm,i,iter);

			if iterCase == 0;     
				V(i) = Viter(2);
				break;						% convergence acheived
			elseif iterCase == 1; break;	% no convergence before maxIter
			elseif iterCase == 2; 
				Viter(1) = Viter(2); 		% continue iterating
			end

		end 		% end iteration loop
	end				% end time loop
end

