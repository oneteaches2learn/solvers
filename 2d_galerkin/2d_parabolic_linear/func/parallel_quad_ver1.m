function int = parallel_quad_ver1(dom,f)
% inputs
%	dom	Domain2d object
%   f	function_handle

	% get edge midpoints
	edges  = dom.getMeshEdges;
	midpts = dom.getMeshEdgeMidpoints(edges);

	% compute F on edge midpoints 
	F = f(midpts(:,1),midpts(:,2));

	% store F as sparse, symmetric, weighted adjacency matrix
	F = sparse([edges(:,1);edges(:,2)], ...
			[edges(:,2);edges(:,1)], ...
			[F(:);F(:)]);

	% loop on elements
	elemNodes = dom.Mesh.Elements';
	elemNodes(:,4) = elemNodes(:,1);
	elemInt = zeros(dom.nElem,1);
	parfor i = 1:dom.nElem

		% loop on edges
		tot = 0;
		for j = 1:3
			tot = tot + F(elemNodes(i,j),elemNodes(i,j+1));
		end
		
		% compute three pt quadrature on element
		elemInt(i) = dom.elemAreas(i) * tot / 3;

	end

	int = sum(elemInt);

end
