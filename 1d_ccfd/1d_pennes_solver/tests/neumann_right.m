function result = neumann_right

	% True Solution
	syms x;
	uTrue = -x * (x - 2) + 2;

	% Domain
	L = 2;

	% Solver Parameters
	k = 2;
	c = 1;
	uStar = 0;

	% Boundary Conditions
	leftBC = 0;
	leftType = 'D';
	rightBC = 0;
	rightType = 'N';

	% Assemble inputs
	domain     = {L,0};
	parameters = {k,c,uStar};
	boundary   = {leftBC,leftType,rightBC,rightType};

	% Test
	result = mms_test(uTrue,domain,parameters,boundary);

end

