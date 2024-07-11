function [BC_L,BCtype_L,BC_R,BCtype_R] = set_subdomain_BCs(boundary,punctures)

	leftBC    = boundary{1};
	leftType  = boundary{2};
	rightBC   = boundary{3};
	rightType = boundary{4};

	subdom_num = length(punctures) + 1;

	% Set BCs for each subdomain
	BC_L{1}     = leftBC;
	BCtype_L{1} = leftType;
	for i = 2:subdom_num
		BC_R{i-1} = punctures{i-1}.BC_L;
		BCtype_R{i-1} = punctures{i-1}.BCtype_L;
		BC_L{i}   = punctures{i-1}.BC_R;
		BCtype_L{i} = punctures{i-1}.BCtype_R;
	end
	BC_R{subdom_num} = rightBC;
	BCtype_R{subdom_num} = rightType;

end
