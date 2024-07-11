function boundary = manufacture_BCs(u,q,a,b,boundary)

	[leftBC,leftType,rightBC,rightType] = unpack_mms_boundary(boundary);

	% manufacture boundary conditions
	if leftType == 'D'
		leftBC = @(t)(u(a,t));
	elseif leftType == 'N'
		leftBC = @(t)(-q(a,t));
	elseif leftType == 'R'
		alpha   = 1;
		sigma   = @(t)(u(a,t) + 1/alpha * q(a,t));
		leftBC  = {alpha,sigma};
	end
	if rightType == 'D'
		rightBC = @(t)(u(b,t));
	elseif rightType == 'N'
		rightBC = @(t)(q(b,t));
	elseif rightType == 'R'
		alpha   = 1;
		sigma   = @(t)(u(a,t) - 1/alpha * q(b,t));
		rightBC  = {alpha,sigma};
	end

	boundary = {leftBC,leftType,rightBC,rightType};
end
