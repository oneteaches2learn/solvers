classdef NewtonGalerkinSolver2d

    properties
    end

    methods
        function self = NewtonGalerkinSolver2d()
            ...
        end

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
					[bCond,t,U] = self.checkVariables(bCond);
                    
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
							1/2 * edgeLength * bCond(coords(edge(1),1),coords(edge(1),2),t,U(edge(1)));
						E(edge(2),edge(2)) = E(edge(2),edge(2)) + ...
							1/2 * edgeLength * bCond(coords(edge(2),1),coords(edge(2),2),t,U(edge(2)));
					end
				end
			end
		end

		function E = computeNonlinearRobinContribution(self)

			%{
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
					[bCond,t,U] = self.checkVariables(bCond);
					[alpha,t,U] = self.checkVariables(alpha);
                    
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
							1/2 * edgeLength * bCond(coords(edge(1),1),coords(edge(1),2),t,sum(U(edge))/2);
						E(edge(2),edge(2)) = E(edge(2),edge(2)) + ...
							1/2 * edgeLength * bCond(coords(edge(2),1),coords(edge(2),2),t,sum(U(edge))/2);

					end
				end
			end
			%}
			nNodes = self.domain.mesh.nNodes;
			E = sparse(zeros(nNodes,nNodes));
			%full(max(max(E)))
		end

    end
end

    