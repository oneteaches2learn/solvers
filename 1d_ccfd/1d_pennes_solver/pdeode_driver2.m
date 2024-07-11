function pdeode_driver2
% notes
%	(1) I want to carefully explain this code snippet:
%
%		% create body/blood temp function
%		c_B   = @(v)(C_B * w_v(v));
%		R     = @(v,C_B)(resolvent(v,dt,C_B,vLower,vUpper));
%		S     = @(u,dx)(weightedAverage(u,dx,w_u));
%		uStar = @(u,v,v_old,dx)(R(dt*mu + dt*S(u,dx) + v_old,c_B(v)));
%
%	c_B is the energy exchange coefficient for the ODE in my coupled system.
%	However, c_B = c_B(v), is given by
%
%		c_B(v) = C_B * w_v(v)
%
%	where C_B is a constant and w_v is a weight function. This makes c_B a
%	source of nonlinearity in the system. Therefore, the value of c_B cannot be
%	set until inside the inner loop that resolves the nonlinearity through
%	iteration. Hence, the best we can do is pass c_B as a function_handle, and
%	then determine the value of c_B in the innermost loop.
%
%	However, c_B is actually used to set a specific resolvent R. The way the
%	resolvent function was designed, the user is to choose dt, C_B, vLower,
%	vUpper in advance, and then create an anonymous function R = R(v) using the
%	syntax
%
%		R = @(v)(resolvent(v,dt,c_B,vLower,vUpper);
%
%	But since the specific resolvent depends on the choice of c_B, and
%	c_B cannot be chosen until the inner iteration loop, we also cannot pick a
%	specific R until the inner iteration loop either. So, we instead set R =
%	@(v,C_B)(...), implying that R is a function of v but is parameterized by
%	C_B.
%
%	Now also, the weightedAverage function was designed to generate anonymous
%	functions of the form S = S(u) by choosing a specific dx and anonymous
%	weight function w, then using the syntax
%
%		S = @(u)(weightedAverage(u,dx,w);
%
%	However, dx will not be set until we are inside the loop, so again, we
%	create instead S = @(u,dx)(...), thereby allowing to parameterize S by dx,
%	once dx has been created. 
%
%	Finally, all of these things are combined to give
%
%		uStar = @(u,v,v_old,dx)(R(dt*mu + dt*S(u,dx) + v_old,c_B(v)));
%
%	The original nonlinear solver expects blood/body temperature to be given by
%	an anonymous function uStar. So, we are essentially just providing the
%	solver what it expects, i.e. an anonymous uStar. It's just that uStar will
%	eventually be evaluated using Uiter(1), Viter(1), and V(i-1), i.e. the
%	previous iteration's U and V, and also V from the previous time step. 
%-----------------------------------------------------------------------------%

	syms x; syms t; syms u; syms v;
	% USER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Domain
	a = 0;
	b = 1;
	M = 100;

	% Solver Parameters
	k     = 1;
	por   = 1;
	f     = 0;

	% cutoff functions
	uLower = 0;
	uUpper = 1; 
	vLower = 0;
	vUpper = 1;
	w_u = @(u)(weightFunction(u,uLower,uUpper));
	w_v = @(v)(weightFunction(v,vLower,vUpper));

	% Tissue energy exchange function
	C = 1;
	%c_T = @(u,v)(C + 0*u + 0*v);
	c_T = @(u,v)(energyExchangeFunction(u,v,w_u,w_v,C));

	% Time Stepping
	dt  = 100;
	nt  = 20;
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
	leftBC = {@()(alpha),@()(v_0)};
	leftType = 'R';
	rightBC = {@()(alpha),@()(0.0)};
	rightType = 'R';

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
	c_B   = @(v)(C_B * w_v(v));
	%c_B   = @(v)(C_B);
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

	[U,xc,tGrid,V] = pennesCCFD1d_sequentialODE(domain,parameters,boundary,time,iterParm,odeParm);


	% Post Processing
	dx = (xc(2) - xc(1)) * ones(length(xc),1);
	for i = 1:length(tGrid), S_unweighted(i) = weightedAverage(U(i,:)',dx,w_Itotal); end
	for i = 1:length(tGrid), S_weighted(i) = weightedAverage(U(i,:)',dx,w_Iweighted); end

	subplot(2,2,1)
	surfplotterCCFD(U,{xc,tGrid},dt,16); 

	subplot(2,2,2)
	plot(tGrid,V)
	title("Blood/Body Temperature")
	xlabel("time")
	ylabel("temperature")

	subplot(2,2,3)
	plot(tGrid,S_unweighted)
	title("Average Hand Temp (whole hand)")
	xlabel("time")
	ylabel("temperature")
	
	subplot(2,2,4)
	plot(tGrid,S_weighted)
	title("Weighted Average Hand Temp")
	xlabel("time")
	ylabel("temperature")
end








