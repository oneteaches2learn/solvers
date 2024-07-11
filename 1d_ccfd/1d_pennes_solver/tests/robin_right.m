function result = robin_right

	% True Solution
	syms x;
	uTrue = -x * (x - 2) + 1;

	% Domain
	L = 2;

	% Solver Parameters
	k = 2;
	c = 0;
	uStar = 0;

	% Boundary Conditions
	alpha = 2;
	leftBC = 1;
	leftType = 'D';
	rightBC = [alpha,0];
	rightType = 'R';

	% Test Parameters
	base          = 5;
	startPower    = 1;
	trials        = 3;
	expectedOrder = 2;
	tolerance     = 0.1;
	displayPlot   = 1;
	displayData   = 1;

	% Assemble inputs
	domain     = {L,0};
	parameters = {k,c,uStar};
	boundary   = {leftBC,leftType,rightBC,rightType};
	testParams = {base,startPower,trials,expectedOrder,tolerance,displayPlot,displayData};

	% Test
	result = mms_test(uTrue,domain,parameters,boundary);

end

