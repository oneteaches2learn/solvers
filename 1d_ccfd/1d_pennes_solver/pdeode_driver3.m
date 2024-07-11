function S = pdeode_driver3

	syms x; syms t; syms u; syms v;
	% USER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Domain
	a = 0;
	b = 1;
	M = 1000;
	h = (b - a) / M;

	% Solver Parameters
	k     = 0.5;
	por   = 10^3;
	f     = 0;

	% cutoff functions
	uLower = 0;
	uUpper = 1; 
	vLower = 0;
	vUpper = 1;
	w_u = @(u)(weightFunction(u,uLower,uUpper));
	w_v = @(v)(weightFunction(v,vLower,vUpper));

	% Tissue energy exchange function
	C = 10^3;
	c_T = @(u,v)(C + 0*u + 0*v);
	%c_T = @(u,v)(energyExchangeFunction(u,v,w_u,w_v,C));

	% Time Stepping
	dt  = 0.1;
	nt  = 1;
	u_0 = 1;

	% Iteration Parameters
	atol    = 1e-6;
	maxIter = 100;
	verb    = 1;

	% ODE Parameters
	v_0 = 1;
	mu  = 0;
	C_B = 1;
	rLower = NaN;
	rUpper = NaN;
	
	% Boundary Conditions
	alpha = 1.0;
	%leftBC = {@()(alpha),@()(v_0)};
	%leftType = 'R';
	%rightBC = {@()(alpha),@()(0.0)};
	%rightType = 'R';
	leftBC = @()(v_0);
	leftType = 'D';
	rightBC = @()(0.0);
	rightType = 'D';

	% Weight for integral average
	weighted_average = 0; 
	w_Iweighted = w_u;
	w_Itotal    = @(u)(1.0 + 0*u);

	% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% ensure input functions have all the same variables
	k     = checkSymbolicFunction(k,x,t,u,v);
	por   = checkSymbolicFunction(por,x,t,u,v);
	f     = checkSymbolicFunction(f,x,t,u,v);
	u_0   = checkSymbolicFunction(u_0,x);

	% convert symbolic functions to function_handles
	k   = matlabFunction(k);
	por = matlabFunction(por);
	f   = matlabFunction(f);
	u_0 = matlabFunction(u_0);

	% create body/blood temp function
	%c_B   = @(v)(C_B * w_v(v));
	c_B   = @(v)(C_B);
	R     = @(v,C_B)(resolvent(v,dt,C_B,rLower,rUpper));
	if weighted_average == 1
		w_I = w_Iweighted;
	else
		w_I = w_Itotal;
	end
	S     = @(u,dx)(weightedAverage(u,dx,w_I));	
	uStar = @(u,v,v_old,dx)(R(dt*mu + dt*c_B(v)*S(u,dx) + v_old,c_B(v)));


	% Assemble Inputs
	domain     = {a,b,M};
	parameters = {k,c_T,uStar,f,por};
	boundary   = {leftBC,leftType,rightBC,rightType};
	time       = {dt,nt,u_0};
	iterParm   = {atol,maxIter,verb};
	odeParm    = {v_0,mu,c_B,vLower,vUpper};

	[U,xc,tGrid,V,S] = pennesCCFD1d_sequentialODE(domain,parameters,boundary,time,iterParm,odeParm);

	full(S)
	eigenvalues = eigs(S,M);
	eig_max = max(eigenvalues)
	eig_min = min(eigenvalues)
	cond = eig_max / eig_min

end








