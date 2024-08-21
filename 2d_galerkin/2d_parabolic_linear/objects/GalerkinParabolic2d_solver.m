classdef GalerkinParabolic2d_solver

	properties
		domain
		time
		coefficients
		uInit
		f
		solution
	end

	properties (Hidden)
		tensors
		vectors
		timestep
		t = 0
	end

	properties (Hidden,Dependent)
		isFirstTimeStep
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

			% initialize problem
			FreeNodes = self.domain.freeNodes;
			self = self.initializeProblem;

			for timestep = 1:self.time.M_t

				% initialize timestep
				self.timestep = timestep;
				self = self.initializeTimestep;

				% assemble problem
				self  = self.assembleTensors;
				self  = self.assembleVectors;
				self  = self.assembleBCs;
				[S,b] = self.finalAssembly;

				% solve and store solution
				v = sparse(self.domain.nNodes,1);
				v(FreeNodes) = S(FreeNodes,FreeNodes) \ b(FreeNodes);
				self.solution(:,self.timestep) = v + self.vectors.U_D;

				% break at equilibrium
				if self.equilibrium == 1, break; end

			end


			%self.solution = U;
			self = self.cleanup;

		end

		function self = initializeProblem(self)

			% Record first timestep
			coords = self.domain.Mesh.Nodes';
			self.solution = zeros(self.domain.nNodes,self.time.N_t);
			self.solution(:,1) = self.uInit(coords(:,1),coords(:,2));

			% Initialize tensors
			self = self.initializeTensors;
			self = self.initializeVectors;

		end

		function self = initializeTensors(self)

			% create fields for tensor storage
			tensors.A   = [];
			tensors.M_p = [];
			tensors.E   = [];
			tensors.M_p_prevTime = self.assembleMassMatrix(self.coefficients.p);

			% check which tensors are time-varying
			tensors.timeVarying.A   = Coefficients.isTimeVarying(self.coefficients.k);
			tensors.timeVarying.M_p = Coefficients.isTimeVarying(self.coefficients.p);

			% update tensors property
			self.tensors = tensors;

		end

		function self = initializeVectors(self)

			% create fields for vector storage
			self.vectors.b_vol = [];
			self.vectors.U_D   = [];
			self.vectors.b_neu = []; 
			self.vectors.b_rob = []; 
			self.vectors.U_prevTime = []; 

			% check which vectors are time-varying
			self.vectors.timeVarying.b_vol = Coefficients.isTimeVarying(self.f);

		end

		function self = initializeTimestep(self)

			% store solution at previous timestep
			self.vectors.U_prevTime = self.solution(:,self.timestep);

			% update time stepping
			self.t = self.timestep * self.time.dt;
			self.timestep = self.timestep + 1;

		end

		function self = assembleTensors(self)

			% if first timestep, create tensors
			if self.isFirstTimeStep == 1
				self.tensors.A   = self.assembleStiffnessMatrix;
				self.tensors.M_p = self.assembleMassMatrix(self.coefficients.p);

			% else, update tensors as needed
			else

				% update A
				if self.tensors.timeVarying.A == 1
					self.tensors.A = self.assembleStiffnessMatrix;
				end

				% update M_p
				if self.tensors.timeVarying.M_p == 1
					self.tensors.M_p_prevTime = self.tensors.M_p;
					cof = self.coefficients.p;
					self.tensors.M_p = self.assembleMassMatrix(cof);
				end

			end


		end

		function self = assembleVectors(self)

			% if first timestep, create vectors
			if self.isFirstTimeStep == 1
				self.vectors.b_vol = self.computeVolumeForces;

			% else, update vectors as needed
			else

				% update b_vol
				if self.vectors.timeVarying.b_vol == 1
					self.vectors.b_vol = self.computeVolumeForces;
				end

			end

		end

		function self = assembleBCs(self)

			self.vectors.U_D   = self.computeDirichletBCs;
			self.vectors.b_neu = self.computeNeumannBCs;
			[self.tensors.E,self.vectors.b_rob] = self.computeRobinBCs;

		end

		function [S,b] = finalAssembly(self)

			% NOTE: placeholder function. Actually handled by specific subclasses.
			...

		end

		function val = get.isFirstTimeStep(self)

			val = (self.t == self.time.dt);

		end

		function A = assembleStiffnessMatrix(self)

			% store variables
			k = self.coefficients.k;

			% check coefficient variables
			if Coefficients.isTimeVarying(k) == 0
				k = @(x1,x2,t)(k(x1,x2));
			end

			% unpack variables
			nNodes    = size(self.domain.Mesh.Nodes,2);
			nElem3    = size(self.domain.Mesh.Elements,2);
			coords    = self.domain.Mesh.Nodes';
			elements3 = self.domain.Mesh.Elements';

			% initialize storage
			A = sparse(nNodes,nNodes);

			% compute k on nodes
			for j = 1:nNodes
				K(j) = k(coords(j,1),coords(j,2),self.t);
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

		function B = assembleMassMatrix(self,c)

			% check coefficient variables
			if Coefficients.isTimeVarying(c) == 0
				c = @(x1,x2,t)(c(x1,x2));
			end
				
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
					self.domain.elemAreas(j) * [2,1,1;1,2,1;1,1,2] / 12;
			end

			% compute r on nodes
			for j = 1:nNodes
				C(j) = c(coords(j,1),coords(j,2),self.t);
			end

			% scale mass matrix by r values
			B = B.*C';

		end

		function M = stima3(self,vertices)

			d = size(vertices,2);
			G = [ones(1,d+1);vertices'] \ [zeros(1,d);eye(d)];
			M = det([ones(1,d+1);vertices']) * G * G' / prod(1:d);

		end

		function b = computeVolumeForces(self)

			% store variables
			f = self.f;

			% check coefficient variables
			if Coefficients.isTimeVarying(f) == 0
				f = @(x1,x2,t)(f(x1,x2));
			end

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
					self.domain.elemAreas(j) * ...
					f(sum(elementCoord(:,1))/3,sum(elementCoord(:,2))/3,self.t) / 3;
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
										coords(bNodes_i,1),coords(bNodes_i,2),self.t);

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
							edgeLength * bCond(edgeMidPt(1),edgeMidPt(2),self.t) / 2;
					end
				end
			end

		end

		function [E,b_rob] = computeRobinBCs(self)

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
							edgeLength * alpha(edgeMidPt(1),edgeMidPt(2),self.t) * ...
							u_R(edgeMidPt(1),edgeMidPt(2),self.t) / 2;

						% compute E matrix
						E(edge(1),edge(1)) = E(edge(1),edge(1)) + ...
							1/2 * edgeLength * alpha(coords(edge(1),:),self.t);
						E(edge(2),edge(2)) = E(edge(2),edge(2)) + ...
							1/2 * edgeLength * alpha(coords(edge(2),:),self.t);

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

		function result = equilibrium(self)

			% default value is zero
			result = 0;

			try % <~~~ if equilibrium is not set, then the code block below
				%	   produces an error. Using the "try" statement allows to
				%	   skip the following block of code in this case. 

			% check if solver should break
			if strcmp(self.time.equilibrium.atEq,"break") == 1
				
				% get error between subsequent timesteps
				arg1 = self.solution(:,self.timestep);
				arg2 = self.solution(:,self.timestep-1);
				errVec = arg1 - arg2;
				err = self.domain.L2norm_threePointQuadrature_nodal(errVec);

				% if error < tolerance, break
				if err < self.time.equilibrium.tolerance
					result = 1;
				end

				% issue `no equilibrium' warning
				if self.timestep == self.time.N_t
					warn = " WARNING, Trial terminated without reaching equilibrium, ";
					fprintf(warn);
				end
				
			end
			end

		end

		function self = cleanup(self)

			self.tensors = [];
			self.vectors = [];

			% update stored parameters
			self.time.N_t = self.timestep;
			self.time.M_t = self.time.N_t - 1;
			self.time.T   = self.time.M_t * self.time.dt;
			self.solution = self.solution(:,1:self.time.N_t);
			self.time.equilibrium.t_eq = self.time.M_t * self.time.dt;

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

		function h = plotPatch(self,timestep)

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

			h = gcf;

		end

		function p = plotMaxSolutionValue(self)

			% store variables
			maxVal = self.getMaxSolutionValue;
			timeGrid = self.time.getTimeGrid;

			% plot result
			p = plot(timeGrid,maxVal);

		end

		function p = plotAverageSolutionValue(self)

			% store variables
			avgVal = self.getAverageSolutionValue;
			timeGrid = self.time.getTimeGrid;

			% plot result
			p = plot(timeGrid,avgVal);

		end

		function animate(self)

			% store variables
			N_t = self.time.N_t;

			% capture zLim bounds
			u_min = min(min(self.solution));
			u_max = max(max(self.solution));
			if u_max == u_min, u_max = u_min + 1; end

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
				self.plotPatch(i);
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

		function result = getAverageSolutionValue(self,timestep)

			% capture average value at each timestep
			avgVal = self.getSolutionIntegral;
			avgVal = avgVal / self.domain.domainArea;

			% if average value at specified timestep is desired
			if nargin == 2
					
				% return average value at specified timestep					
				result = avgVal(timestep);

			% else return vector of average values at each timestep
			else
				result = avgVal;
			end

		end

		function int = getSolutionIntegral(self,NameValueArgs)

			arguments
				self
				NameValueArgs.timestep double = 0
				NameValueArgs.quadOrder double = 1
			end

			% if no timestep passed, loop over all time steps
			if NameValueArgs.timestep == 0
				time_start = 1;
				time_stop = self.time.N_t;

			% else compute at a specific timestep
			else
				time_start = NameValueArgs.timestep;
				time_stop = NameValueArgs.timestep;
			end

			% call quadrature
			for i = time_start:time_stop

				int(i) = self.domain.nodalQuadrature(self.solution(:,i));
			
			end
	
		end

		%{
		% DEPRECATED! ---------------------------------------------------------%
		function int = getsolutionintegral(self,NameValueArgs)

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
			
				% loop over timesteps
				%int = zeros(1,time_stop - time_start + 1);
				parfor n = time_start:time_stop

					% store numerical solution
					U_n = self.solution(:,n);

					% loop over elements
					ind = n - time_start + 1;
					int(n) = 0;
					for j = 1:self.domain.nElem
						elemInd = self.domain.Mesh.Elements(:,j);
						int(n) = int(n) + (sum(U_n(elemInd)) / 3) * ...
							self.domain.elemAreas(j);
					end
					%temp(ind) = tot;

				end
				int = 0;

			% else run slower algorithm for higher order quadrature
			else
				int = zeros(1,time_stop - time_start + 1);
				for n = time_start:time_stop

					% store numerical solution
					U_n = self.solution(:,n);

					% loop over elements
					nElem = size(self.domain.Mesh.Elements,2);
					ind = n - time_start + 1;
					for i = 1:nElem

						% interpolate u_h locally 
						locNodes = self.domain.Mesh.Elements(:,i); 
						locCoord = self.domain.Mesh.Nodes(:,locNodes)';
						uLoc = scatteredInterpolant(locCoord,U_n(locNodes)); 

						% compute local error function on quad points
						Qp = quadtriangle(NameValueArgs.quadOrder,'Domain',locCoord); 
						uSol_Qp = uLoc(Qp.Points(:,1),Qp.Points(:,2)); 

						% quadrature on local error function
						int(ind) = int(ind) + dot(Qp.Weights, uSol_Qp); 

					end
				end
			end
	
		end
		%----------------------------------------------------------------------%
		%}

		%{
		% DEPRECATED! Use L2 inner products in Domain2d instead----------------%
		function IP = L2_IP(self,arg1,arg2,NameValueArgs)

			arguments
				self
				arg1
				arg2
				NameValueArgs.quadOrder double = 1
			end

			% if first order quadrature, run much faster algorithm
			if NameValueArgs.quadOrder == 0

				% compute product of arguments
				arg = arg1 .* arg2;
		
				% loop over elements
				IP = 0;
				for j = 1:self.domain.nElem
					elemInd = self.domain.Mesh.Elements(:,j);
					IP = IP + (sum(arg(elemInd)) / 3) * ...
						self.domain.elemAreas(j);
				end

			% else run slower algorithm for higher order quadrature
			else

				% loop over elements
				IP = 0;
				for i = 1:self.domain.nElem

					% interpolate arguments locally 
					locNodes = self.domain.Mesh.Elements(:,i); 
					locCoord = self.domain.Mesh.Nodes(:,locNodes)';
					arg1_interp = scatteredInterpolant(locCoord,arg1(locNodes)); 
					arg2_interp = scatteredInterpolant(locCoord,arg2(locNodes)); 

					% compute arguments on local quad points
					Qp = quadtriangle(NameValueArgs.quadOrder,'Domain',locCoord); 
					arg1_Qp = arg1_interp(Qp.Points(:,1),Qp.Points(:,2)); 
					arg2_Qp = arg2_interp(Qp.Points(:,1),Qp.Points(:,2)); 

					% quadrature on product of arguments
					IP = IP + dot(Qp.Weights,arg1_Qp .* arg2_Qp); 

				end
			end

		end
		%---------------------------------------------------------------------%
		%}

		function [t,timestep] = getConvergenceTime(self,tol)

			% default tolerance = 10^-6
			if nargin == 1, tol = 10^-6; end

			for n = 1:self.time.M_t

				% store current and next timestep
				arg1 = self.solution(:,n);		
				arg2 = self.solution(:,n+1);		

				% compute error between subsequent timesteps
				err_n = arg1 - arg2;

				% compute L2 norm of error
				int = self.domain.L2norm_piecewiseLinear(err_n);

				% check against tolerance
				if int < tol
					timestep = n-1;
					t = self.time.dt * timestep;
					return	
				end

			end

			% if loop completes, then the desired tolerance is never reached
			t = NaN;
			timestep = NaN;
			fprintf('WARNING: Tolerance not reached during time-series\n');
		end

	end
end
