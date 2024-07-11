classdef GalerkinPoisson2d_solver

	properties
		domain
		coefficients
		f
		solution
	end

	methods
		function self = GalerkinPoisson2d_solver(dom,cofs,f)
			
			% store inputs
			self.domain = dom;
			self.coefficients = cofs;
			self.f = f;
			
			% calculate solution
			self = self.solve;

		end

		function self = solve(self)

			% load variables
			coords = self.domain.Mesh.Nodes';
			FreeNodes = self.domain.freeNodes;
			r = self.coefficients.r;
			k = self.coefficients.k;

			% assemble tensors
			A = self.assembleStiffnessMatrix;
			B = self.assembleMassMatrix;

			% update RHS and BCs
			b_vol = self.computeVolumeForces();
			U_D = self.computeDirichletBCs;
			b_neu = self.computeNeumannBCs;
			[E,b_rob] = self.computeRobinBCs;

			% final assembly
			v = sparse(self.domain.nNodes,1);
			S = A + B + E;
			b = b_vol - b_neu + b_rob - S * U_D;

			% solve and store solution
			v(FreeNodes) = S(FreeNodes,FreeNodes) \ b(FreeNodes);
			self.solution = v + U_D;

		end

		function A = assembleStiffnessMatrix(self)

			% unpack variables
			nNodes    = size(self.domain.Mesh.Nodes,2);
			nElem3    = size(self.domain.Mesh.Elements,2);
			coords    = self.domain.Mesh.Nodes';
			elements3 = self.domain.Mesh.Elements';
			k = self.coefficients.k;

			% initialize storage
			A = sparse(nNodes,nNodes);

			% compute k on nodes
			for j = 1:nNodes
				K(j) = k(coords(j,1),coords(j,2));
			end

			% assemble stiffness matrix
			for j = 1:nElem3
				elementInd   = elements3(j,:);
				elementCoord = coords(elementInd,:);
				K_j = (1/3) * sum(K(elementInd));
				A(elementInd,elementInd) = A(elementInd,elementInd) + ...
				 	K_j * self.stima3(elementCoord);
			end

		end

		function B = assembleMassMatrix(self)

			% unpack variables
			nNodes    = size(self.domain.Mesh.Nodes,2);
			nElem3    = size(self.domain.Mesh.Elements,2);
			coords 	  = self.domain.Mesh.Nodes';
			elements3 = self.domain.Mesh.Elements';
			r = self.coefficients.r;

			% initialize storage
			B = sparse(nNodes,nNodes);
			
			% assemble base mass matrix
			for j = 1:nElem3
				elementInd = elements3(j,:);
				elementCoord = coords(elementInd,:);
				B(elementInd,elementInd) = B(elementInd,elementInd) + ...
					det([1,1,1;elementCoord']) * [2,1,1;1,2,1;1,1,2] / 24;
			end

			% compute r on nodes
			for j = 1:nNodes
				R(j) = r(coords(j,1),coords(j,2));
			end

			% scale mass matrix by r values
			B = B.*R';

		end

		function M = stima3(self,vertices)

			d = size(vertices,2);
			G = [ones(1,d+1);vertices'] \ [zeros(1,d);eye(d)];
			M = det([ones(1,d+1);vertices']) * G * G' / prod(1:d);

		end

		function b = computeVolumeForces(self)

			% unpack variables
			nNodes    = size(self.domain.Mesh.Nodes,2);
			nElem3    = size(self.domain.Mesh.Elements,2);
			coords    = self.domain.Mesh.Nodes';
			elements3 = self.domain.Mesh.Elements';

			% initialize storage
			b = sparse(nNodes,1);

			% compute volume forces
			for j = 1:nElem3
				elementInd    = elements3(j,:);
				elementCoord  = coords(elementInd,:);
				b(elementInd) = b(elementInd) + ...
					det([1,1,1; elementCoord']) * ...
					self.f(sum(elementCoord(:,1))/3,sum(elementCoord(:,2))/3) / 6;
			end
		end

		function U_D = computeDirichletBCs(self)

			% unpack variables
			dom    = self.domain;
			nNodes = self.domain.nNodes;
			coords = self.domain.Mesh.Nodes';

			% initialize storage
			U_D = sparse(nNodes,1);

			% compute Dirichlet boundary conditions
			for i = 1:dom.NumEdges
				
				if dom.edges(i).boundaryType == 'D'

					% store nodes on i-th edge of domain
					bNodes_i = dom.edges(i).nodes;

					% compute value of U on nodes of i-th edge
					U_D(bNodes_i) = dom.edges(i).boundaryCondition(...
										coords(bNodes_i,1),coords(bNodes_i,2));

				end
			end

		end

		function b_neu = computeNeumannBCs(self)

			% unpack variables
			dom    = self.domain;
			nNodes = self.domain.nNodes;
			coords = self.domain.Mesh.Nodes';

			% initialize storage
			b_neu = sparse(nNodes,1);

			% compute boundary conditions
			for i = 1:self.domain.NumEdges
				
				% compute Neumann condition
				if dom.edges(i).boundaryType == 'N'

					% store nodes on i-th edge of domain
					bNodes_i = dom.edges(i).nodes;
					bCond = dom.edges(i).boundaryCondition;

					% loop over segments of i-th edge
					for j = 1:length(bNodes_i)-1
						edge = [bNodes_i(j) bNodes_i(j+1)];
						edgeMidPt = sum(coords(edge,:)/2);
						edgeLength = norm(coords(edge(1),:) - coords(edge(2),:));
						b_neu(edge) = b_neu(edge) + ...
							edgeLength * bCond(edgeMidPt(1),edgeMidPt(2)) / 2;
					end
				end
			end

		end


		function [E,b_rob] = computeRobinBCs(self,BCnodes)

			% unpack variables
			dom    = self.domain;
			nNodes = self.domain.nNodes;
			coords = self.domain.Mesh.Nodes';

			% initialize storage
			b_rob = sparse(nNodes,1);
			E = sparse(nNodes,nNodes);

			% compute boundary conditions
			for i = 1:self.domain.NumEdges
				
				% compute Dirichlet condition
				if dom.edges(i).boundaryType == 'R'
					
					bCond = dom.edges(i).boundaryCondition;

					% unpack functions
					alpha = bCond{1};
					u_R = bCond{2};

					% store nodes on i-th edge of domain
					bNodes_i = dom.edges(i).nodes;

					% loop over segments of i-th edge
					for j = 1:length(bNodes_i)-1

						% get edge data
						edge = [bNodes_i(j) bNodes_i(j+1)];
						edgeMidPt = sum(coords(edge,:)/2);
						edgeLength = norm(coords(edge(1),:) - coords(edge(2),:));

						% compute RHS vector
						b_rob(edge) = b_rob(edge) + ...
							edgeLength * alpha(edgeMidPt(1),edgeMidPt(2)) * ...
							u_R(edgeMidPt(1),edgeMidPt(2)) / 2;

						% compute E matrix
						E(edge(1),edge(1)) = E(edge(1),edge(1)) + ...
							1/2 * edgeLength * alpha(coords(edge(1),:));
						E(edge(2),edge(2)) = E(edge(2),edge(2)) + ...
							1/2 * edgeLength * alpha(coords(edge(2),:));

					end
				end
			end
		end

		function bcNodes = getBCnodes(self)

			fegm = self.domain;
			gd   = self.domain.geometryMatrix;

			% loop over objects used to define domain
			midPts = [];
			for i = 1:size(gd,2) 

				% get object nodes
				nodes = [];
				for j = 1:4
					nodes = [nodes; gd(j+2,i) gd((j+4)+2,i)];
				end

				% add midpoints to list
				for j = 0:3
					midPts = [midPts; (nodes(j+1,:) + nodes(mod(j+1,4)+1,:)) / 2;];
				end
			end
			
			% collect boundary node pairs
			bcNodes = [];
			for i = 1:length(midPts)
				edgeID = fegm.nearestEdge(midPts(i,:));
				nodes = fegm.Mesh.findNodes('region','Edge',edgeID);
				for j = 1:length(nodes)-1
					bcNodes = [bcNodes; nodes(j) nodes(j+1)];
				end
			end
		end

		function plot(self)

			% store domain information
			coordinates = self.domain.Mesh.Nodes';
			elements3 = self.domain.Mesh.Elements';
			elements4 = [];

			% get solution at final time step
			U = full(self.solution);

			% plot data
			trisurf(elements3,coordinates(:,1),coordinates(:,2),U', ...
				'facecolor','interp')
			hold on
			trisurf(elements4,coordinates(:,1),coordinates(:,2),U', ...
				'facecolor','interp')
			hold off

			% format plot
			view(10,40);

		end

	end
end
