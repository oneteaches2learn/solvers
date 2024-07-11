function A_h = computeMassMatrix(por,dx)
%CONSTRUCTMASSMATRIX(POR,DX) constructs diagonal mass matrix using porosity.

	% instantiate storage
	A_h = sparse(length(por),length(por));

	% add entries on diagonal
	for i = 1:length(por)
		A_h(i,i) = por(i) * dx(i);
	end
end
