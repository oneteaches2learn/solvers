function result = mms_driver_nonlinear
% MMS_DRIVER_NONLINEAR drives the function mms_test_nonlinear to perform a set
%	 of mms tests for the pennesCCFD1dTimeNonlinear solver
%
% Author: Tyler Fara					Date: April 2, 2024
%-----------------------------------------------------------------------------%
% Inputs
%	None
%
% Outputs
% 	result		binary, 1 = all tests passed, 0 = at least one test failed
%-----------------------------------------------------------------------------%
% Notes
%	(1) Treat MMS_DRIVER_NONLINEAR like a template: Duplicate it, rename the
%	duplicates, encapsulate them with specific projects worked on, etc. 
%
%	(2) The user sets the settings under USER PARAMETERS below. Treat the current
%	settings like examples, and change the settings as desired, following these
%	examples. The function uTrue is the desired true solution. The parameters
%	a, b represent the domain endpoints. The parameters k, c, uStar, por are
%	(respectively) conductivity, Helmholtz coefficient, blood temperature, and
%	porosity. boundaryTypes 'DNR' represent the types of boundary conditions
%	that should be tested
%
%	(3) uTrue and coefficients can be typed in exactly as they appear on some
%	piece of paper. That is, if c = x^2 + 0.5t + u, then type in c = x^2 +
%	0.5*t + u. If uStar = 1, then type in uStar = 1. After being passed to
%	mms_test_nonlinear, any function F will be converted to a corresponding
%	symfun F(x,t,u). This allows to manipulate the functions using MATLAB's
%	symbolic function engine in order to manufacture the necessary boundary
%	conditions, initial condition, and RHS. Then, these symbolic inputs will be
%	converted to function_handles and passed to the solver; the solver uses
%	function_handle coefficients because evaluating function_handles is much
%	faster than evaluating symbolic functions. 
%
%	(4) As of April 3, 2024, nonlinearities in every coefficient *except* k can
%	be handled. So, for now, pennesCCFD1dTimeNonlinear is only a *semilinear*
%	solver. At some point, I might want to try to extend it to a quasilinear
%	solver (i.e. one that can handle nonlinearities in k). But for my current
%	purpose, I only need a semilinear solver. To the user: k = k(x,t) is fine.
%	k = k(x,t,u) may not work. 
%-----------------------------------------------------------------------------%

	% USER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% True Solution
	syms x; % do not change
	syms t; % do not change
	syms u; % do not change
	uTrue = sin(pi*x) * exp(-t) + t;

	% Domain
	a = 1;
	b = 2;

	% Solver Parameters
	k     = 1*x + t;
	c     = x^2 + 0.5 + t + u;
	uStar = x - 2 + t + u;
	por   = x + 1 + t + u;

	% Iteration Parameters
	iterParm = 'default';

	% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Boundary Conditions
	boundaryTypes = 'DNR';

	tic	
	wait = waitbar(0,'Please wait...');
	for i = 1:3, for j = 1:3
		leftType  = boundaryTypes(i);
		rightType = boundaryTypes(j);
	
		% Assemble inputs
		domain     = {a,b,0};
		parameters = {k,c,uStar,por};
		boundary   = {0,leftType,0,rightType};

		% Test
		trial_num = (i-1)*3+j;
		test_results(trial_num) = mms_test_nonlinear(uTrue,domain,parameters,boundary,iterParm,'default');
		
		waitbar(trial_num / 9,wait,'Please wait...')

		fprintf('Boundary Type %c%c, result: %i\n',leftType,rightType,test_results(trial_num));
	end, end
	toc
	
	close(wait);
	result = min(test_results);

	if result == 0
		test_results
	end

end
