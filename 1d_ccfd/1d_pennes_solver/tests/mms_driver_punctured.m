function result = mms_driver_punctured
% MMS_DRIVER_PUNCTURED drives the function mms_test_punctured to perform a set
% 	of mms tests for the pennesCCFD1dtime_punctured solver
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
%	(1) Treat MMS_DRIVER_PUNCTURED like a template: Duplicate it, rename the
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
%	(3) All functions should be symbolic functions. After being passed to
%	mms_test_time, these symbolic functions will be manipulated in order to
%	produce boundary conditions, initial conditions, and a right hand side that
%	manufactures the desired solution. These manipulations require symbolic
%	functions. 
%
%	(4) The input punctures should be a cell array of Puncture objects.
%	Puncture objects are instantiated as 
%	
%		punc = Puncture(L,R,L_BC,R_BC,L_type,R_type);
%
%	where L and R are the left and right endpoints, L_BC, R_BC are the left and
%	right boundary conditions, and L_type and R_type are 'D', 'N', or 'R' and
%	represent the boundary condition type. Note that mms_test_punctured will
%	manufacture L_BC and R_BC based on uTrue. So these inputs will be
%	overwritten, hence the @()(0.0) defaults. 
%-----------------------------------------------------------------------------%

	% USER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% True Solution
	syms x; % do not change
	syms t; % do not change
	uTrue(x,t) = sin(pi*x) * exp(-t) + t;

	% Domain
	a = 0;
	b = 1;

	% Solver Parameters
	k     = 10*x + 1;
	c     = x^2 + 0.5;
	uStar = x - 2;
	por   = x + 1;

	% Punctures
	punc1 = Puncture(0.2,0.3,@()(0.0),@()(0.0),'D','D');
	punc2 = Puncture(0.7,0.8,@()(0.0),@()(0.0),'D','D');
	punctures = {punc1,punc2};

	% Settings for MMS test
	mmsParms = 'default';


	% BLACK BOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
		test_results(trial_num) = mms_test_punctured(uTrue,domain,parameters,boundary,mmsParms,punctures);
		
		waitbar(trial_num / 9,wait,'Please wait...')
	end, end
	toc
	
	close(wait);
	result = min(test_results);

	if result == 0
		test_results
	end
end

