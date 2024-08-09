function int = parallel_quad_ver2(dom,f,N_t)
% inputs
%	dom	Domain2d object
%   f	function_handle
%	N_t integer, number of timesteps

	% get edge midpoints
	edges  = dom.getMeshEdges;
	midpts = dom.getMeshEdgeMidpoints(edges);

	% get element data
	elemNodes = dom.Mesh.Elements';
	elemNodes(:,4) = elemNodes(:,1);

	% loop on timesteps
	int = zeros(1,N_t);
	parfor n = 1:N_t

		% compute F on edge midpoints 
		F = f(midpts(:,1),midpts(:,2));

		% store F as sparse, symmetric, weighted adjacency matrix
		F = sparse([edges(:,1);edges(:,2)], ...
				[edges(:,2);edges(:,1)], ...
				[F(:);F(:)]);

		% loop on elements
		for i = 1:dom.nElem

			% loop on edges
			tot = 0;
			for j = 1:3
				tot = tot + F(elemNodes(i,j),elemNodes(i,j+1));
			end
			
			% compute three pt quadrature on element
			int(n) = int(n) + dom.elemAreas(i) * tot/3;

		end

	end

end

