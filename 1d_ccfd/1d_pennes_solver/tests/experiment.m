function result = experiment

	% True Solution
	syms x;
	uTrue = sin(pi * x) - 1;

	% Domain
	L = 1.5;

	% Solver Parameters
	k = 1;
	c = 10;
	uStar = -3;

	% Boundary Conditions
	leftBC = 0;
	leftType = 'N';
	rightBC = 0;
	rightType = 'N';

	% Assemble inputs
	domain     = {L,0};
	parameters = {k,c,uStar};
	boundary   = {leftBC,leftType,rightBC,rightType};

	% Test
	result = mms_test(uTrue,domain,parameters,boundary);

end

