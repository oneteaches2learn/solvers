function [leftBC,leftType,rightBC,rightType] = unpack_mms_boundary(boundary)
	leftBC    = boundary{1};
	leftType  = boundary{2};
	rightBC   = boundary{3};
	rightType = boundary{4};
end
