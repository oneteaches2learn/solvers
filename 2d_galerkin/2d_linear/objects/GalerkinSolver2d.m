classdef GalerkinSolver2d

	properties
		domain
		coefficients
		f
		solution
	end

	properties (Hidden)
		tensors
		vectors
	end

	methods
        % CONSTRUCTOR
		function self = GalerkinSolver2d(dom,auxfun)

			if nargin == 2
				
				% store inputs
				self.domain = dom;
				self.coefficients = auxfun.cofs;
				self.f = auxfun.f;

			end
			
			% calculate solution
			% ... 

		end


        % ASSEMBLY FUNCTIONS
		function A = assembleStiffnessMatrix(self)

			% store variables
			k = self.coefficients.k;

			% check coefficient variables
			if Coefficients.isTimeVarying(k) == 0
				k = @(x1,x2,t)(k(x1,x2));
                t = 0;
            else
                t = self.t;
			end

			% unpack variables
			nNodes    = self.domain.mesh.nNodes;
			nElem3    = self.domain.mesh.nElems;
			coords    = self.domain.mesh.nodes;
			elements3 = self.domain.mesh.elements;

			% initialize storage
			A = sparse(nNodes,nNodes);

			% compute k on nodes
			K = k(coords(:,1),coords(:,2),t);
			if size(K,1) > 1
				K_avg = sum(K(elements3),2) / 3;
			else
				K_avg = K * ones(self.domain.mesh.nElems,1);
			end

			% assemble stiffness matrix
			for j = self.domain.mesh.effectiveElems;
				elementInd   = elements3(j,:);
				elementCoord = coords(elementInd,:);
				A(elementInd,elementInd) = A(elementInd,elementInd) + ...
					K_avg(j) * self.stima3(elementCoord);

			end
		end

		function B = assembleMassMatrix(self,c)

			% unpack variables
			nNodes    = self.domain.mesh.nNodes;
			nElem3    = self.domain.mesh.nElems;
			coords 	  = self.domain.mesh.nodes;
			elements3 = self.domain.mesh.elements;

			% check coefficient variables
			[c,t,U] = self.checkVariables(c);
				
			% initialize storage
			B = sparse(nNodes,nNodes);
			
			% assemble base mass matrix
			for j = self.domain.mesh.effectiveElems
				elementInd = elements3(j,:);
				elementCoord = coords(elementInd,:);
				B(elementInd,elementInd) = B(elementInd,elementInd) + ...
					self.domain.mesh.areas(j) * [2,1,1;1,2,1;1,1,2] / 12;
			end

			% compute c on nodes
			C(:) = c(coords(:,1),coords(:,2),t,U);

			% scale mass matrix by c values
			B = B.*C';

		end

		function M = stima3(self,vertices)

			d = size(vertices,2);
			G = [ones(1,d+1);vertices'] \ [zeros(1,d);eye(d)];
			M = det([ones(1,d+1);vertices']) * G * G' / prod(1:d);

		end

		function b = computeVolumeForces(self,f)

			% unpack variables
			nNodes    = self.domain.mesh.nNodes;
			nElem3    = self.domain.mesh.nElems;
			coords    = self.domain.mesh.nodes;
			elements3 = self.domain.mesh.elements;

			% store function
			if nargin == 1,
				f = self.f;
			end

			% check coefficient variables
			[f,t,U] = self.checkVariables(f);

			% initialize storage
			b = sparse(nNodes,1);

			% compute volume forces
			for j = self.domain.mesh.effectiveElems
				elementInd    = elements3(j,:);
				elementCoord  = coords(elementInd,:);
				b(elementInd) = b(elementInd) + ...
					self.domain.mesh.areas(j) * ...
					f( sum(elementCoord(:,1))/3, sum(elementCoord(:,2))/3, t, sum(U(elementInd))/3 ) / 3;
			end
		end

		function [f,t,U] = checkVariables(self,f)

			% unpack variables
			nNodes    = self.domain.mesh.nNodes;

			% check coefficient variables
			if ~Coefficients.isNonlinear(f) && ~Coefficients.isTimeVarying(f)
				f = @(x1,x2,t,u)(f(x1,x2));
				t = 0;
				U = zeros(nNodes,1);
			elseif ~Coefficients.isNonlinear(f) && Coefficients.isTimeVarying(f)
				f = @(x1,x2,t,u)(f(x1,x2,t));
				t = self.t;
				U = zeros(nNodes,1);
			elseif Coefficients.isNonlinear(f) && ~Coefficients.isTimeVarying(f)
				f = @(x1,x2,t,u)(f(x1,x2,u));
				t = 0;
				U = self.U;
			else
				t = self.t;
				U = self.U;
			end
		end

		function self = assembleBCs(self)

			self = self.computePeriodicBCs;
			self.vectors.U_D   = self.computeDirichletBCs;
			self.vectors.b_neu = self.computeNeumannBCs;
			[self.tensors.M_rob,self.vectors.b_rob] = self.computeRobinBCs;
			[self.tensors.M_rob,self.tensors.M_dyn,self.vectors.b_dyn] = self.computeDynamicBCs;

		end


        % BOUNDARY CONDITIONS
		function U_D = computeDirichletBCs(self)

			% unpack variables
			dom    = self.domain;
			nNodes = self.domain.mesh.nNodes;
			coords = self.domain.mesh.nodes;

			% initialize storage
			U_D = sparse(nNodes,1);

			% compute Dirichlet boundary conditions
			for i = 1:dom.boundary.nEdges
				
				if dom.boundary.edges(i).boundaryType == 'D'

					% store nodes on i-th edge of domain
					bNodes_i = dom.boundary.edges(i).nodes;

                    % check if boundary condition is time-varying
                    if Coefficients.isTimeVarying(dom.boundary.edges(i).boundaryCondition) == 0
                        dom.boundary.edges(i).boundaryCondition = @(x1,x2,t)(dom.boundary.edges(i).boundaryCondition(x1,x2));
                        t = 0;
                    else
                        t = self.t;
                    end

					% compute value of U on nodes of i-th edge
					U_D(bNodes_i) = dom.boundary.edges(i).boundaryCondition(...
										coords(bNodes_i,1),coords(bNodes_i,2),t);

				end
			end
		end

		function b_neu = computeNeumannBCs(self)

			% unpack variables
			dom    = self.domain;
			nNodes = self.domain.mesh.nNodes;
			coords = self.domain.mesh.nodes;

			% initialize storage
			b_neu = sparse(nNodes,1);

			% compute boundary conditions
			for i = 1:self.domain.boundary.nEdges
				
				% compute Neumann condition
				if dom.boundary.edges(i).boundaryType == 'N'

					% store nodes on i-th edge of domain
					bNodes_i = dom.boundary.edges(i).nodes;
					bCond = dom.boundary.edges(i).boundaryCondition;

					%{
                    % OLD: check if boundary condition is time-varying
                    if Coefficients.isTimeVarying(bCond) == 0
                        bCond = @(x1,x2,t)(bCond(x1,x2));
                        t = 0;
                    else
                        t = self.t;
                    end
					%}

					% NEW: check boundary condition variables
					[bCond,t,U] = checkVariables(self,bCond);

					% loop over segments of i-th edge
					for j = 1:length(bNodes_i)-1
						edge = [bNodes_i(j) bNodes_i(j+1)];
						edgeMidPt = sum(coords(edge,:)/2);
						edgeLength = norm(coords(edge(1),:) - coords(edge(2),:));
						b_neu(edge) = b_neu(edge) + ...
							edgeLength * bCond(edgeMidPt(1),edgeMidPt(2),t,sum(U(edge))/2) / 2;
					end
				end
			end
		end

		function [E,b] = computeRobinBCs(self)

			% unpack variables
			dom    = self.domain;
			nNodes = self.domain.mesh.nNodes;
			coords = self.domain.mesh.nodes;

			% initialize storage
			b = sparse(nNodes,1);
			E = sparse(nNodes,nNodes);

			% compute boundary conditions
			for i = 1:self.domain.boundary.nEdges
				
				% compute Dirichlet condition
				if dom.boundary.edges(i).boundaryType == 'R'
					
					bCond = dom.boundary.edges(i).boundaryCondition;

					% unpack functions
					alpha = bCond{1};
					u_BC  = bCond{2};

                    % check if alpha is time-varying
                    if Coefficients.isTimeVarying(alpha) == 0
                        alpha = @(x1,x2,t)(alpha(x1,x2));
                        t = 0;
                    else
                        t = self.t;
                    end

                    % check if u_R is time-varying
                    if Coefficients.isTimeVarying(u_BC) == 0
                        u_BC = @(x1,x2,t)(u_BC(x1,x2));
                        t = 0;
                    else
                        t = self.t;
                    end
                    
					% store nodes on i-th edge of domain
					bNodes_i = dom.boundary.edges(i).nodes;

					% loop over segments of i-th edge
					for j = 1:length(bNodes_i)-1

						% get edge data
						edge = [bNodes_i(j) bNodes_i(j+1)];
						edgeMidPt = sum(coords(edge,:)/2);
						edgeLength = norm(coords(edge(1),:) - coords(edge(2),:));

						% compute RHS vector
						b(edge) = b(edge) + ...
							edgeLength * alpha(edgeMidPt(1),edgeMidPt(2),t) * ...
							u_BC(edgeMidPt(1),edgeMidPt(2),t) / 2;

						% compute E matrix
						E(edge(1),edge(1)) = E(edge(1),edge(1)) + ...
							1/2 * edgeLength * alpha(coords(edge(1),:),t);
						E(edge(2),edge(2)) = E(edge(2),edge(2)) + ...
							1/2 * edgeLength * alpha(coords(edge(2),:),t);

					end
				end
			end
		end

		function [E,F,g] = computeDynamicBCs(self)

			% unpack variables
			dom    = self.domain;
			nNodes = self.domain.mesh.nNodes;
			coords = self.domain.mesh.nodes;

			% initialize storage
			E = sparse(nNodes,nNodes);
			F = sparse(nNodes,nNodes);
			g = sparse(nNodes,1);

			% compute boundary conditions
			for i = 1:self.domain.boundary.nEdges
				
				% compute Dirichlet condition
				if dom.boundary.edges(i).boundaryType == 'T'
					
					bCond = dom.boundary.edges(i).boundaryCondition;

					% unpack functions
					alpha = bCond{1};
					beta  = bCond{2};
					gamma = bCond{3};
					u_BC  = bCond{4};

                    % check if alpha is time-varying
                    if Coefficients.isTimeVarying(alpha) == 0
                        alpha = @(x1,x2,t)(alpha(x1,x2));
                        t = 0;
                    else
                        t = self.t;
                    end

                    % check if beta is time-varying
                    if Coefficients.isTimeVarying(beta) == 0
                        beta = @(x1,x2,t)(beta(x1,x2));
                        t = 0;
                    else
                        t = self.t;
                    end

                    % check if alpha is time-varying
                    if Coefficients.isTimeVarying(gamma) == 0
                        gamma = @(x1,x2,t)(gamma(x1,x2));
                        t = 0;
                    else
                        t = self.t;
                    end

                    % check if u_R is time-varying
                    if Coefficients.isTimeVarying(u_BC) == 0
                        u_BC = @(x1,x2,t)(u_BC(x1,x2));
                        t = 0;
                    else
                        t = self.t;
                    end
                    
					% store nodes on i-th edge of domain
					bNodes_i = dom.boundary.edges(i).nodes;

					% loop over segments of i-th edge
					for j = 1:length(bNodes_i)-1

						% get edge data
						edge = [bNodes_i(j) bNodes_i(j+1)];
						edgeMidPt = sum(coords(edge,:)/2);
						edgeLength = norm(coords(edge(1),:) - coords(edge(2),:));

						% E matrix corresponds to u term
						E(edge(1),edge(1)) = E(edge(1),edge(1)) + ...
							1/2 * edgeLength * alpha(coords(edge(1),:),t);
						E(edge(2),edge(2)) = E(edge(2),edge(2)) + ...
							1/2 * edgeLength * alpha(coords(edge(2),:),t);

						% F matrix corresponds to du/dt term
						F(edge(1),edge(1)) = F(edge(1),edge(1)) + ...
							1/2 * edgeLength * beta(coords(edge(1),:),t);
						F(edge(2),edge(2)) = F(edge(2),edge(2)) + ...
							1/2 * edgeLength * beta(coords(edge(2),:),t);

						% compute RHS vector
						g(edge) = g(edge) + ...
							edgeLength * gamma(edgeMidPt(1),edgeMidPt(2),t) * ...
							u_BC(edgeMidPt(1),edgeMidPt(2),t) / 2;

					end
				end
			end
		end

		function self = computePeriodicBCs(self)

			% unpack variables
			dom    = self.domain;
			nNodes = self.domain.mesh.nNodes;
			coords = self.domain.mesh.nodes;
			P = self.domain.boundary.P_nodes;

			% get periodic source term
			b_per = sparse(nNodes,1);
			b_vol = self.vectors.b_vol;
			b_per(P.free.edge) = b_vol(P.replica.edge);
			b_per(P.free.corner) = sum(b_vol(P.replica.corner));
			self.vectors.b_per = b_per;

			% get periodic matrix
			F = sparse(nNodes,nNodes);
			F_rows = sparse(nNodes,nNodes);
			F_cols = sparse(nNodes,nNodes);
			F_diag = sparse(nNodes,nNodes);

			% copy replica rows of tensor to corresponding free rows
			F_rows(P.free.edge,:) = self.tensors.A(P.replica.edge,:);
			if length(P.free.corner) > 0
				F_rows(P.free.corner,:) = sum(self.tensors.A(P.replica.corner,:),1);
			end
			
			% copy replica cols of tensor to corresponding free cols
			F_cols(:,P.free.edge) = self.tensors.A(:,P.replica.edge);
			if length(P.free.corner) > 0
				F_cols(:,P.free.corner) = sum(self.tensors.A(:,P.replica.corner),2);
			end
			
			% copy replica diag entries of tensor to corresponding free diag entries
			F_diag(P.free.edge,:) = F_cols(P.replica.edge,:);
			if length(P.free.corner) > 0
				F_diag(P.free.corner,:) = sum(F_cols(P.replica.corner,:),1);
			end

			% assemble final result
			F = F_rows + F_cols + F_diag;
			self.tensors.P = F;

		end

		function bcNodes = getBCnodes(self)

			fegm = self.domain;
			gd   = self.domain.geometryMatrix;

			% loop over objects used to define domain
			midPts = [];
			for i = 1:size(gd,2) 

				% get object nodes
				nodes = []
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
				nodes = fegm.mesh.Mesh.findNodes('region','Edge',edgeID);
				for j = 1:length(nodes)-1
					bcNodes = [bcNodes; nodes(j) nodes(j+1)];
				end
			end
		end


        % UTILITY FUNCTIONS
		function self = cleanup(self)

            % clear out tensors and vectors
			self.tensors = [];
			self.vectors = [];

			% add NaN for unused nodes
			self.solution(self.domain.mesh.unusedNodes,:) = NaN;

			% copy solution to periodic replica nodes
			P = self.domain.boundary.P_nodes;
			self.solution(P.replica.edge,:) = self.solution(P.free.edge,:);
			self.solution(P.replica.corner,:) = ...
							repmat(self.solution(P.free.corner,:),3,1);					
			%self.solution(P.replica.corner,:) = self.solution(P.free.corner,:);

            % ensure solution is full
            self.solution = full(self.solution);

		end


		% PLOTTING FUNCTIONS
		function h = plot(self,timestep)

			% for parabolic problems, plot final time point unless otherwise specified
			if isa(self,'GalerkinSolver2d_parabolic') && nargin < 2, 
                timestep = self.domain.time.N_t; 
            elseif isa(self,'GalerkinSolver2d_elliptic')
                timestep = 1;
            end

			% store domain information
			coordinates = self.domain.mesh.nodes;
			elements3 = self.domain.mesh.elements;
			elements4 = [];

			% get solution at desired time step
			U = full(self.solution(:,timestep));

			% plot data
			trisurf(elements3,coordinates(:,1),coordinates(:,2),U', ...
				'facecolor','interp')
			hold on
			trisurf(elements4,coordinates(:,1),coordinates(:,2),U', ...
				'facecolor','interp')
			hold off

			% format plot
			view(10,40);
			%title('Solution of the Problem')

			h = gcf;
		end

		function h = plotPatch(self,timestep)

			% for parabolic problems, plot final time point unless otherwise specified
			if isa(self,'GalerkinSolver2d_parabolic') && nargin < 2, 
                timestep = self.domain.time.N_t; 
            elseif isa(self,'GalerkinSolver2d_elliptic')
                timestep = 1;
            end

			% store domain information
			coordinates = self.domain.mesh.nodes;
			elements3 = self.domain.mesh.elements;

			% get solution at desired time step
			u_min = min(min(self.solution));
			u_max = max(max(self.solution));
			U = full(self.solution(:,timestep));

			% plot data
			patch('Faces',elements3, ...
				  'Vertices',coordinates, ...
				  'FaceVertexCData',U, ...
				  'FaceColor','interp');
			clim([u_min u_max])
			colorbar

			h = gcf;

		end

		function plot_yline(self,NameValuePairs)
		% for now, assumes only 1 yline

			arguments
				self
				NameValuePairs.plotInclusions = 'off'
				NameValuePairs.timestep = 0;
			end

			% store variables
			dl = self.domain.boundary.dl;
			mesh = self.domain.mesh.Mesh;

			% get yline nodes
			segIDs1 = dl.segIDs_yLines;
			segIDs2 = dl.segIDs_yLines_inc;
			nodes1 = findNodes(mesh,"region","Edge",segIDs1);
			nodes2 = findNodes(mesh,"region","Edge",segIDs2);
			nodes = [nodes1, nodes2];

			% get x- and solution vals; sort by x vals
			x = self.domain.mesh.nodes(nodes,1);
			sol = self.solution(nodes,:);
			if NameValuePairs.timestep == 0
				sol = sol(:,end);
			else
				sol = sol(:,NameValuePairs.timestep);
			end
			val = sortrows([x,sol]);

			% get yline height
			y = dl.mat(dlObj.yLowRow,segIDs1(1));

			% plot
			plot(val(:,1),val(:,2),"LineWidth",6);
			
			% set yLim
			uMax = max(max(self.solution));
			uMin = min(min(self.solution));
			ylim([uMin uMax*1.1]);

			% get x-vals at edges of inclusions
			if strcmp(NameValuePairs.plotInclusions,'on')
				if ~isempty(segIDs2)
					for i = 1:length(segIDs2)
						nodes_temp = findNodes(mesh,"region","Edge",segIDs2(i));
						incBounds(2*i-1) = self.domain.mesh.nodes(nodes_temp(1),1);
						incBounds(2*i) = self.domain.mesh.nodes(nodes_temp(end),1);
					end
				end
				
				% plot inclusion bounds
				xline(incBounds,'--');
			end

		end


		% VTK EXPORT FUNCTIONS
		function tria2vtk(self,filename,dir,timesteps)

			% prepare directory
			if isempty(dir), dir = cd; end
			if ~strcmpi(dir(end),'/'); dir(end+1)='/'; end

			% store data
			nodes = self.domain.mesh.nodes;
			elems = self.domain.mesh.elements;
			data = self.solution;

			if nargin == 3
				vtk.tria2vtk(filename,nodes,elems,data(:,end),dir);
			else
				for i = 1:length(timesteps)
					if timesteps(i) > self.domain.time.N_t
						timesteps(i) = self.domain.time.N_t;
					end
					filename_i = sprintf('%s_%i',filename,i);
					vtk.tria2vtk( ...
						filename_i,nodes,elems,data(:,timesteps(i)),dir);
				end
			end

		end

		function tria2vtk_mesh(self,filename,dir)

			% prepare directory
			if isempty(dir), dir = cd; end
			if ~strcmpi(dir(end),'/'); dir(end+1)='/'; end

			% store data
			nodes = self.domain.mesh.nodes;
			elems = self.domain.mesh.elements;

			filename = sprintf('%s_mesh',filename);
			vtk.tria2vtk_mesh(filename,nodes,elems,dir);

		end

		function tria2vtk_surf(self,filename,dir,timesteps)

			% prepare directory
			if isempty(dir), dir = cd; end
			if ~strcmpi(dir(end),'/'); dir(end+1)='/'; end

			% store data
			nodes = self.domain.mesh.nodes;
			elems = self.domain.mesh.elements;
			data = self.solution;

			if nargin == 3
				vtk.tria2vtk_surf(filename,nodes,elems,data(:,end),dir);
			else
				for i = 1:length(timesteps)
					if timesteps(i) > self.domain.time.N_t
						timesteps(i) = self.domain.time.N_t;
					end
					filename_i = sprintf('%s_surf_%i',filename,i);
					vtk.tria2vtk_surf(...
						filename_i,nodes,elems,data(:,timesteps(i)),dir);
				end
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
			avgVal = avgVal / self.domain.getDomainArea;

			% if average value at specified timestep is desired
			if nargin == 2
					
				% return average value at specified timestep					
				result = avgVal(timestep);

			% else return vector of average values at each timestep
			else
				result = avgVal;
			end
		end

		function int = getSolutionIntegral(self,timestep)

            % reshape solution for integration
            input(:,1,:) = self.solution(:,:);

            % if timestep is specified, integrate at that timestep
            if nargin == 2
			    int = self.domain.nodalQuadrature(input(:,1,timestep));

            % else integrate across all timesteps
            else
                int = self.domain.nodalQuadrature(input);
            end

		end

		function result = getMaxSolutionValue_omegaEps(self,timestep)

			% get solution on omegaEps
			sol = self.solution_omegaEps;

			% capture max value at each timestep
			maxVal = max(sol);

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

		function result = getMinSolutionValue_omegaEps(self,timestep)

			% get solution on omegaEps
			sol = self.solution_omegaEps;

			% capture min value at each timestep
			minVal = min(sol);

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

		function result = getAverageSolutionValue_omegaEps(self,timestep)

			% capture average value at each timestep
			avgVal = self.getSolutionIntegral_omegaEps;
			avgVal = avgVal / (self.domain.getDomainArea * self.domain.inclusion.volumeFraction);

			% if average value at specified timestep is desired
			if nargin == 2
					
				% return average value at specified timestep					
				result = avgVal(timestep);

			% else return vector of average values at each timestep
			else
				result = avgVal;
			end
		end

		function int = getSolutionIntegral_omegaEps(self,timestep)

			% get solution on omegaEps
			sol = self.solution_omegaEps;

            % reshape solution for integration
            input(:,1,:) = sol(:,:);

            % if timestep is specified, integrate at that timestep
            if nargin == 2
			    int = self.domain.nodalQuadrature(input(:,1,timestep));

            % else integrate across all timesteps
            else
                int = self.domain.nodalQuadrature(input);
            end

		end

		function sol = solution_omegaEps(self)

			keep_nodes = [];
			self.domain.boundary.dl.faceIDs_Omega_eps
			for i = self.domain.boundary.dl.faceIDs_Omega_eps
				keep_nodes = [keep_nodes, self.domain.mesh.Mesh.findNodes('region','Face',i)];
			end

			discard_nodes = setdiff(1:self.domain.mesh.nNodes,keep_nodes);

			sol = self.solution;
			sol(discard_nodes,:) = NaN;

		end



    end
end