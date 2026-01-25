classdef NewtonGalerkinSolver1d

    properties
    end

	properties (Hidden)
		DJ  % Jacobian matrix
		J   % Residual vector
	end

    methods
        function self = NewtonGalerkinSolver1d()
            ...
        end

		function E = computeNonlinearNeumannContribution(self)

			% unpack variables
			dom    = self.domain;
			nNodes = dom.mesh.nNodes;
			coords = dom.mesh.nodes;   % assumed nNodes x 1 or nNodes x dim with dim=1

			% initialize storage
			E = sparse(nNodes, nNodes);

			% loop over boundary entities
			for i = 1:dom.boundary.nEdges

				% only treat Neumann boundaries
				if dom.boundary.edges(i).boundaryType == 'N'

					% nonlinear Neumann term: d(g)/d(u) (or similar)
					bCond = dom.boundary.edges(i).boundaryCondition_ddu;

					% check if time-varying etc.
					[bCond, t, U, V] = self.checkVariables(bCond);

					% nodes on this Neumann boundary (in 1D: just points)
					bNodes_i = dom.boundary.edges(i).nodes;

					% loop over boundary nodes on this edge
					for k = 1:length(bNodes_i)

						node = bNodes_i(k);
						x    = coords(node, 1);  % 1D coordinate

						% contribution is purely local and diagonal at this boundary node
						E(node, node) = E(node, node) + ...
							bCond(x, t, U(node), V);

					end
				end
			end
		end

		% note: dummy implementation; real feature not yet implemented in 1d
		function E = computeNonlinearRobinContribution(self)

			E = sparse(self.domain.mesh.nNodes, self.domain.mesh.nNodes);

		end



		%{
		% note: 2d version kept for reference
        function E = computeNonlinearNeumannContribution(self)

			% unpack variables
			dom    = self.domain;
			nNodes = self.domain.mesh.nNodes;
			coords = self.domain.mesh.nodes;

			% initialize storage
			E = sparse(nNodes,nNodes);

			% compute boundary conditions
			for i = 1:self.domain.boundary.nEdges
				
				% compute Dirichlet condition
				if dom.boundary.edges(i).boundaryType == 'N'
					
					bCond = dom.boundary.edges(i).boundaryCondition_ddu;

                    % check if alpha is time-varying
					[bCond,t,U,V] = self.checkVariables(bCond);
                    
					% store nodes on i-th edge of domain
					bNodes_i = dom.boundary.edges(i).nodes;

					% loop over segments of i-th edge
					for j = 1:length(bNodes_i)-1

						% get edge data
						edge = [bNodes_i(j) bNodes_i(j+1)];
						edgeMidPt = sum(coords(edge,:)/2);
						edgeLength = norm(coords(edge(1),:) - coords(edge(2),:));

						% compute E matrix
						E(edge(1),edge(1)) = E(edge(1),edge(1)) + ...
							1/2 * edgeLength * bCond(coords(edge(1),1),coords(edge(1),2),t,U(edge(1)),V);
						E(edge(2),edge(2)) = E(edge(2),edge(2)) + ...
							1/2 * edgeLength * bCond(coords(edge(2),1),coords(edge(2),2),t,U(edge(2)),V);
					end
				end
			end
		end
		%}

		%{
		% note: 2d version kept for reference
		function E = computeNonlinearRobinContribution(self)

			% unpack variables
			dom    = self.domain;
			nNodes = self.domain.mesh.nNodes;
			coords = self.domain.mesh.nodes;

			% initialize storage
			E = sparse(nNodes,nNodes);

			% compute boundary conditions
			for i = 1:self.domain.boundary.nEdges
				
				% compute Dirichlet condition
				if dom.boundary.edges(i).boundaryType == 'R'
					
					alpha = dom.boundary.edges(i).boundaryCondition{1};
					bCond = dom.boundary.edges(i).boundaryCondition_ddu;

                    % check if alpha is time-varying
					[bCond,t,U,V] = self.checkVariables(bCond);
					[alpha,t,U,V] = self.checkVariables(alpha);
                    
					% store nodes on i-th edge of domain
					bNodes_i = dom.boundary.edges(i).nodes;

					% loop over segments of i-th edge
					for j = 1:length(bNodes_i)-1

						% get edge data
						edge = [bNodes_i(j) bNodes_i(j+1)];
						edgeMidPt = sum(coords(edge,:)/2);
						edgeLength = norm(coords(edge(1),:) - coords(edge(2),:));

						% compute E matrix
						E(edge(1),edge(1)) = E(edge(1),edge(1)) + ...
							1/2 * edgeLength * bCond(coords(edge(1),1),coords(edge(1),2),t,sum(U(edge))/2,V);
						E(edge(2),edge(2)) = E(edge(2),edge(2)) + ...
							1/2 * edgeLength * bCond(coords(edge(2),1),coords(edge(2),2),t,sum(U(edge))/2,V);

					end
				end
			end
		end
		%}

    end
end

    