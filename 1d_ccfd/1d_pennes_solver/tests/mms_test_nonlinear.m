function result = mms_test_nonlinear(uTrue,domain,parameters,boundary,iterParm,mmsParms)
% MMS_TEST_TIME(UTRUE,DOMAIN,PARAMETERS,BOUNDARY,MMSPARMS,VARARGIN) performs an
% 	MMS test on the pennesCCFD1dtime solver. 
%
% Author: Tyler Fara						Date: April 2, 2024
%------------------------------------------------------------------------------%
% Inputs
%	uTrue		symbolic function, representing the true solution
%	domain		cell array, {a,b,0}: 
%					a, double, left domain endpoint
%					b, double, right domain endpoint
%					0, integer, # of cells***
%	parameters	cell array, {k,c,uStar,por} 
%					k,     symbolic function, conductivity
%					c,     symbolic function, helmholtz coefficient
%					uStar, symbolic function, blood temperature
%					por,   symbolic function, porosity
%	boundary	cell array, {0,leftType,0,rightType}
%					0,	matlabFunction, left BC***
%					leftType, char, 'D', 'N', or 'R'
%					0,	matlabFunction, right BC***
%					rightType, char, 'D', 'N', or 'R'
%	mmsParms	{base,startPower,trials,expectedOrder,tolerance,displayPlot,displayData}
%					base, integer, base for cell number
%					startPower, integer, starting exponent base^startPower
%					trials, integer, number of increments of exponent
%					expectedOrder, double, expected order of convergence
%					tolerance, double, allowed deviation from expectedOrder
%					displayPlot, binary, 1 = display plot, 0 = not
%					displayData, binary, 1 = display data, 0 = not
%	varargin	not used (used for punctures in mms_test_punctured)	
%	
%	*** designates an input whose value is immaterial because this input
%	 will be overwritten by MMS_TEST_TIME during the mms test.  
%
% Outputs
%	result		binary, 1 = pass, 0 = fail
%
%-----------------------------------------------------------------------------%
% Notes
%	(1) Coefficients input into MMS_TEST_TIME should be symbolic functions.
%	This is because the MMS test will compute partial derivatives u_x and u_t,
%	flux q, u_init, and the desired boundary conditions from these symbolic
%	function inputs. Then it will compute the necessary RHS f used to
%	manufacture the solution. All of this computation requires symbolic
%	manipulation; hence the inputs need to be symbolic functions. By contrast,
%	coefficients input into CCFD1DTIME should be function_handles.  Thus, the
%	MMS_TEST_TIME converts all of these symbolic functions to
%	function_handle obects using the built-in matlabFunction function. 
%
%	(2) the mmsParms argument can be 
%
%			mmsParms = 'default'
%
%		in which case, the default values are:
%			
%			base		  = 10; 
%			startPower    = 1;
%			trials        = 4;
%			expectedOrder = 2; 
%			tolerance     = 0.1;
%			displayPlot   = 0;
%			displayData   = 0;
%
%-----------------------------------------------------------------------------%

	% unpack inputs
	[a,b] = unpack_mms_domain(domain);
	[k,c,uStar,por] = unpack_mms_parameters(parameters);
	[leftBC,leftType,rightBC,rightType] = unpack_mms_boundary(boundary);
	[base,startPower,trials,expectedOrder,tolerance,displayPlot,displayData] = unpack_mmsParms(mmsParms);
	
	% MANUFACTURE SOLUTION
	% define symbolic variables
	syms x;
	syms t;
	syms u;

	% ensure input functions have all the same variables
	k     = checkSymbolicFunction(k,x,t,u);
	c     = checkSymbolicFunction(c,x,t,u);
	uStar = checkSymbolicFunction(uStar,x,t,u);
	por   = checkSymbolicFunction(por,x,t,u);

	% manufacture forcing functions
	uTrue = checkSymbolicFunction(uTrue,x,t,u);
	u_x(x,t,u) = diff(uTrue,x);
	u_t(x,t,u) = diff(uTrue,t);
	q(x,t,u)   = -k * u_x;
	f(x,t,u)  = por * u_t + diff(q) + c * (uTrue - uStar);
	u_0 = uTrue(x,0,u);

	% convert symbolic functions to function_handles
	u   = matlabFunction(uTrue);
	u_x = matlabFunction(u_x);
	u_t = matlabFunction(u_t);
	k   = matlabFunction(k);
	c   = matlabFunction(c);
	uStar = matlabFunction(uStar);
	por = matlabFunction(por);
	q   = matlabFunction(q);
	f   = matlabFunction(f);
	u_0 = matlabFunction(u_0);

	% manufacture BCs
	boundary = manufacture_BCs(u,q,a,b,boundary);

	% assemble inputs
	parameters = {k,c,uStar,f,por};

	% INITIALIZE STORAGE
	xc = {}; tgrid = {}; U = {}; V = {}; Norm = {}; Ratio = {}; Log = {};

	n = trials;
	% X: TEST
	for i = 1:trials,
		nx = 2^(4+(i-1)); nt = 2^(4+2*(i-1)); dx = (b-a)/nx; dt = (b-a)/nt;
		[U{i},xc{i},tGrid{i}] = pennesCCFD1dTimeNonlinear({a,b,nx},parameters,boundary,{1/nt,nt,u_0},iterParm);
		[X,Y] = meshgrid(xc{i},tGrid{i});
		V{i} = u(X,Y);
		Norm{i} = sqrt(sum(dx*dt*(U{i}-V{i}).^2,'all')); % L2 error
	end;

	% X: PRINT RESULTS
	for i = 1:n-1
		Ratio{i} = Norm{i}/Norm{i+1}; Log{i} = log10(Ratio{i})/log10(2);
	end;

	if displayData == 1, 
		boundary
		Norm, Ratio, Log
	end
	if displayPlot == 1, 
		surfplotterCCFD(U{4},{xc{4},tGrid{4}},dt,16); 
	end

	if (abs(Log{trials-1} - expectedOrder) < tolerance), result = 1; else, result = 0; end

end








