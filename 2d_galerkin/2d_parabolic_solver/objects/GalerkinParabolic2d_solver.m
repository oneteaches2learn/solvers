classdef GalerkinParabolic2d_solver

	properties
		domain
		time
		coefficients
		uInit
		f
		solution
	end

	methods
		function self = GalerkinParabolic2d_solver(dom,time,cofs,uInit,f)
			
			% store inputs
			self.domain = dom;
			self.time = time;
			self.coefficients = cofs;
			self.uInit = uInit;
			self.f = f;
			
			% calculate solution
			self = self.solve;

		end

		function self = solve(self)

			% load variables
			coords = self.domain.Mesh.Nodes';
			FreeNodes = self.domain.freeNodes;

			% compute initial condition
			U = zeros(self.domain.nNodes,self.time.N_t);
			U(:,1) = self.uInit(coords(:,1),coords(:,2));

			for n = 1:self.time.M_t

				% current time
				t = n * self.time.dt;

				% assemble problem
				tensors = self.assembleTensors(t);
				vectors = self.assembleVectors(t);
				[S,b]   = self.finalAssembly(tensors,vectors,U(:,n));

				% solve and store solution
				v = sparse(self.domain.nNodes,1);
				v(FreeNodes) = S(FreeNodes,FreeNodes) \ b(FreeNodes);
				U(:,n+1) = v + vectors.U_D;

			end

			self.solution = U;

		end

		function tensors = assembleTensors(self,t)

			% NOTE: placeholder function. Actually handled by specific subclasses.
			...

		end

		function vectors = assembleVectors(self,t)

			% NOTE: placeholder function. Actually handled by specific subclasses.
			...

		end

		function [S,b] = finalAssembly(self,tensors,vectors,U_prev)

			% NOTE: placeholder function. Actually handled by specific subclasses.
			...

		end

		function A = assembleStiffnessMatrix(self,t)

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
				K(j) = k(coords(j,1),coords(j,2),t);
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

		function B = assembleMassMatrix(self,c,t)

			% unpack variables
			nNodes    = size(self.domain.Mesh.Nodes,2);
			nElem3    = size(self.domain.Mesh.Elements,2);
			coords 	  = self.domain.Mesh.Nodes';
			elements3 = self.domain.Mesh.Elements';

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
				C(j) = c(coords(j,1),coords(j,2),t);
			end

			% scale mass matrix by r values
			B = B.*C';

		end

		function M = stima3(self,vertices)

			d = size(vertices,2);
			G = [ones(1,d+1);vertices'] \ [zeros(1,d);eye(d)];
			M = det([ones(1,d+1);vertices']) * G * G' / prod(1:d);

		end

		function b = computeVolumeForces(self,t)

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
					self.f(sum(elementCoord(:,1))/3,sum(elementCoord(:,2))/3,t) / 6;
			end

		end

		function U_D = computeDirichletBCs(self,t)

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
										coords(bNodes_i,1),coords(bNodes_i,2),t);

				end
			end

		end

		function b_neu = computeNeumannBCs(self,t)

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
							edgeLength * bCond(edgeMidPt(1),edgeMidPt(2),t) / 2;
					end
				end
			end

		end


		function [E,b_rob] = computeRobinBCs(self,t)

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
							edgeLength * alpha(edgeMidPt(1),edgeMidPt(2),t) * ...
							u_R(edgeMidPt(1),edgeMidPt(2),t) / 2;

						% compute E matrix
						E(edge(1),edge(1)) = E(edge(1),edge(1)) + ...
							1/2 * edgeLength * alpha(coords(edge(1),:),t);
						E(edge(2),edge(2)) = E(edge(2),edge(2)) + ...
							1/2 * edgeLength * alpha(coords(edge(2),:),t);

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


		% PLOTTING FUNCTIONS
		function plot(self,timestep)

			% plot final time point unless otherwise specified
			if nargin < 2, timestep = self.time.N_t; end

			% store domain information
			coordinates = self.domain.Mesh.Nodes';
			elements3 = self.domain.Mesh.Elements';
			elements4 = [];

			% get solution at desired time step
			U = self.solution(:,timestep);

			% plot data
			trisurf(elements3,coordinates(:,1),coordinates(:,2),U', ...
				'facecolor','interp')
			hold on
			trisurf(elements4,coordinates(:,1),coordinates(:,2),U', ...
				'facecolor','interp')
			hold off

			% format plot
			view(10,40);
			title('Solution of the Problem')
		end

		function plotPatch(self,timestep)

			% plot final time point unless otherwise specified
			if nargin < 2, timestep = self.time.N_t; end

			% store domain information
			coordinates = self.domain.Mesh.Nodes';
			elements3 = self.domain.Mesh.Elements';

			% get solution at desired time step
			u_min = min(min(self.solution));
			u_max = max(max(self.solution));
			U = self.solution(:,timestep);

			% plot data
			patch('Faces',elements3, ...
				  'Vertices',coordinates, ...
				  'FaceVertexCData',U, ...
				  'FaceColor','interp');
			clim([u_min u_max])
			colorbar

		end

		%{
		function plotMaxSolutionValue(self)

			% store variabes
			maxVal = self.getMaxSolutionValue;
			timeGrid = 
		%}

		function animate(self)

			u_min = min(min(self.solution));
			u_max = max(max(self.solution));
			N_t = self.time.N_t;

			% plot first time step
			self.plot(1);
			zlim([u_min u_max])
			pause();

			% animate remaining steps
			for i = 1:N_t
				self.plot(i)
				zlim([u_min u_max])
				pause(1/N_t)
			end

		end

		function animatePatch(self)

			N_t = self.time.N_t;

			% plot first time step
			self.plotPatch(1);
			pause();

			% animate remaining steps
			for i = 1:N_t
				self.plotPatch(i)
				pause(1/N_t)
			end

		end


		% SOLUTION ANALYSIS
		function result = getMaxSolutionValue(self,timestep)

			% capture max value at each timestep
			maxVal = max(self.solution);

			% if max value at specified timestep is desired
			if nargin == 2

				% return max value across all timesteps,
				if strcmp(timestep,'all')
					result = max(maxVal);

				% return max value at given timestep
				else
					result = maxVal(timestep);
				end

			% else return vector of max values at each timestep
			else
				result = maxVal;
			end

		end

		function result = getMinSolutionValue(self,timestep)

			% capture min value at each timestep
			minVal = min(self.solution);

			% if min value at specified timestep is desired
			if nargin == 2

				% return min value across all timesteps,
				if strcmp(timestep,'all')
					result = min(minVal);

				% return min value at given timestep
				else
					result = minVal(timestep);
				end

			% else return vector of min values at each timestep
			else
				result = minVal;
			end

		end

		function int = getSolutionIntegral(self,NameValueArgs)

			arguments
				self
				NameValueArgs.timestep double = 0
				NameValueArgs.quadOrder double = 1
			end

			tic
			% if no timestep passed, loop over all time steps
			if NameValueArgs.timestep == 0
				time_start = 1;
				time_stop = self.time.N_t;

			% else compute at a specific timestep
			else
				time_start = NameValueArgs.timestep;
				time_stop = NameValueArgs.timestep;
			end

			% if first order quadrature, run much faster algorithm
			if NameValueArgs.quadOrder == 1
			
				% compute element areas
				[temp,elemArea] = area(self.domain.Mesh);

				% loop over timesteps
				int = zeros(1,time_stop - time_start + 1);
				for n = time_start:time_stop

					% store numerical solution
					U_n = self.solution(:,n);

					% loop over elements
					for j = 1:self.domain.nElem
						elemInd = self.domain.Mesh.Elements(:,j);
						int(n) = int(n) + (sum(U_n(elemInd)) / 3) * ...
							self.domain.elemAreas(j);
					end

				end

			% else run slower algorithm for higher order quadrature
			else
				int = zeros(1,time_stop - time_start + 1);
				for n = time_start:time_stop

					% store numerical solution
					U_n = self.solution(:,n);

					% loop over elements
					nElem = size(self.domain.Mesh.Elements,2);
					for i = 1:nElem

						% interpolate u_h locally 
						locNodes = self.domain.Mesh.Elements(:,i); 
						locCoord = self.domain.Mesh.Nodes(:,locNodes)';
						uLoc = scatteredInterpolant(locCoord,U_n(locNodes)); 

						% compute local error function on quad points
						Qp = quadtriangle(NameValueArgs.quadOrder,'Domain',locCoord); 
						uSol_Qp = uLoc(Qp.Points(:,1),Qp.Points(:,2)); 

						% quadrature on local error function
						int(n) = int(n) + dot(Qp.Weights, uSol_Qp); 

					end
				end
			end

		end


	end
end
