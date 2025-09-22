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

			A_vec = self.localStiffnessMatrix_vec;
			A = self.sumLocalMatrices(A_vec);

		end

		function M = assembleMassMatrix(self,c)

			M_vec = self.localMassMatrix_vec(c);
			M = self.sumLocalMatrices(M_vec);

		end

		function A_vec = localStiffnessMatrix_vec(self)

			% store variables
			k = self.coefficients.k;

			% unpack variables
			coords    = self.domain.mesh.nodes;
			elements3 = self.domain.mesh.elements(self.domain.mesh.effectiveElems,:);
			centroids = self.domain.mesh.centroids(self.domain.mesh.effectiveElems,:);
			nElems    = size(elements3,1);

			% check coefficient variables
			[k,t,U,V] = self.checkVariables(k);

			% compute k on centroids
			K = k(centroids(:,1),centroids(:,2),t,U,V);
			if length(K) == 1
				K = K * ones(nElems,1);
			end

			% get 2 x nElem arrays of node coordinates by element
			node1 = coords(elements3(:,1),:)';
			node2 = coords(elements3(:,2),:)';
			node3 = coords(elements3(:,3),:)';

			% rotate into 2 x 3 x nElem sized array
			node1 = permute(node1,[1,3,2]);
			node2 = permute(node2,[1,3,2]);
			node3 = permute(node3,[1,3,2]);
			coords_vec = cat(2, node1, node2, node3);

			% pad with 1's to make 3 x 3 x nElem sized array
			pad = ones(1,3,nElems);
			coords_vec = cat(1, pad, coords_vec);

			% compute G_vec and G_vec'
			divisor = [zeros(1,2);eye(2)];
			G_vec = pagemldivide(coords_vec,divisor);
			G_vec_prime = permute(G_vec,[2,1,3]);

			% Compute unscaled A_vec
			A_vec = pagemtimes(G_vec, G_vec_prime);

			% scale by K to complete centroid quadrature
			K_vec = permute(K,[3,2,1]);
			A_vec = pagemtimes(K_vec,A_vec);

		end

		function M_vec = localMassMatrix_vec(self,c);

			% unpack variables
			centroids = self.domain.mesh.centroids(self.domain.mesh.effectiveElems, :);
			nElems    = size(centroids,1);
    		elements3 = self.domain.mesh.elements(self.domain.mesh.effectiveElems, :);

			% check coefficient variables
			[c,t,U,V] = self.checkVariables(c);

			% interpolate U on centroids
			U_centroids = mean(U(elements3), 2);

			% compute c on centroids
			C = c(centroids(:,1),centroids(:,2),t,U_centroids,V);
			if length(C) == 1
				C = C * ones(nElems,1);
			end

			% compute local mass matrices
			M_loc = self.localMassMatrix();
			C_vec = permute(C,[3,2,1]);
			M_vec = pagemtimes(M_loc,C_vec);

		end

		function M = localMassMatrix(self)

			M = [2,1,1; 1,2,1; 1,1,2] / 12;

		end

		function M = sumLocalMatrices(self,M_vec)

			% unpack variables
			nNodes     = self.domain.mesh.nNodes;
			elementInd = self.domain.mesh.elements(self.domain.mesh.effectiveElems,:);

			% scale M_vec by areas to complete centroid quadrature
			areas = self.domain.mesh.areas(self.domain.mesh.effectiveElems);
			areas = permute(areas,[3,2,1]);
			M_vec = pagemtimes(areas,M_vec);

			% Get summation locations
			[I, J] = ndgrid(1:3, 1:3);
			ind_vec = sub2ind([nNodes,nNodes], elementInd(:,I(:)), elementInd(:,J(:)));
			ind_vec = reshape(ind_vec',[],1);

			% sum the values
			% note: flag 'true' below indicates sparse matrix is to be created
			M_vec_flat = M_vec(:);
			M = accumarray(ind_vec,M_vec_flat,[nNodes * nNodes,1],[],[],true);
			M = reshape(M,[nNodes,nNodes]);

		end

		function b = computeVolumeForces(self,f1,f2)

			% unpack variables
			nNodes    = self.domain.mesh.nNodes;
			elements3 = self.domain.mesh.elements(self.domain.mesh.effectiveElems, :);
			nElem3    = size(elements3,1);
			coords    = self.domain.mesh.nodes;

			% store function
			if nargin == 1,
				f1 = self.f;
				f2 = @(x1,x2,t)1; 
			elseif nargin == 2
				f2 = @(x1,x2,t)1;
			end

			% check coefficient variables
			[f1,t,U,V] = self.checkVariables(f1);
			f2 = self.checkVariables(f2);

			% initialize storage
			b = sparse(nNodes,1);

			% compute centroids of elements
			elementCoordsX = reshape(coords(elements3, 1), [nElem3, 3]);
			elementCoordsY = reshape(coords(elements3, 2), [nElem3, 3]);
			centroidsX = sum(elementCoordsX, 2) / 3;
			centroidsY = sum(elementCoordsY, 2) / 3;
			centroidsU = sum(U(elements3), 2) / 3;

			% compute volume forces
			areas = self.domain.mesh.areas(self.domain.mesh.effectiveElems);
			f1_vals = f1(centroidsX, centroidsY, t, centroidsU, V);
			f2_vals = f2(centroidsX, centroidsY, t, centroidsU, V);
			forces = f1_vals .* f2_vals / 3;
			volumeForces = areas .* forces;

			% accumulate forces into the global vector
			b = accumarray(elements3(:), repmat(volumeForces, 3, 1), [nNodes, 1], @sum, 0, true);

		end

		function self = assembleBCs(self)

			%self = self.computePeriodicBCs;
			self.vectors.U_D   = self.computeDirichletBCs;
			self.vectors.b_neu = self.computeNeumannBCs;
			[self.tensors.M_rob,self.vectors.b_rob] = self.computeRobinBCs;
			[self.tensors.M_dyn_u,self.tensors.M_dyn_du,self.vectors.b_dyn] = self.computeDynamicBCs;

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

					% check boundary condition variables
					[bCond,t,U,V] = checkVariables(self,bCond);

					% loop over segments of i-th edge
					for j = 1:length(bNodes_i)-1
						edge = [bNodes_i(j) bNodes_i(j+1)];
						edgeMidPt = sum(coords(edge,:)/2);
						edgeLength = norm(coords(edge(1),:) - coords(edge(2),:));
						b_neu(edge) = b_neu(edge) + ...
							edgeLength * bCond(edgeMidPt(1),edgeMidPt(2),t,sum(U(edge))/2, V) / 2;
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

					% NEW: check boundary condition variables
					[alpha,t,U,V] = checkVariables(self,alpha);
					[u_BC,t,U,V] = checkVariables(self,u_BC);

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
							edgeLength * alpha(edgeMidPt(1),edgeMidPt(2),t,sum(U(edge))/2,V) * ...
							u_BC(edgeMidPt(1),edgeMidPt(2),t,sum(U(edge))/2,V) / 2;

						% compute E matrix
						E(edge(1),edge(1)) = E(edge(1),edge(1)) + ...
							1/2 * edgeLength * alpha(coords(edge(1),1),coords(edge(1),2),t,sum(U(edge))/2);
						E(edge(2),edge(2)) = E(edge(2),edge(2)) + ...
							1/2 * edgeLength * alpha(coords(edge(2),1),coords(edge(2),2),t,sum(U(edge))/2);

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

		%{
		% NOTE: Deprecated. This function turned out to work, but only for
		% problems with a stiffness matrix only, and not with more complicated
		% tensor structures. On June 26, 2025, it was split into
		% computePeriodicCorrectionMatrix and computePeriodicCorrectionVector.

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
		%}

		function [S,b] = periodicCorrection(self,S,b)

			% compute corrections
			S_corr = self.computePeriodicCorrectionMatrix(S);
			b_corr = self.computePeriodicCorrectionVector(b);

			% return results
			S = S + S_corr;
			b = b + b_corr;

		end

		function S_corr = computePeriodicCorrectionMatrix(self,S)

			% unpack variables
			nNodes = self.domain.mesh.nNodes;
			P = self.domain.boundary.P_nodes;

			% initialize storage
			F_rows = sparse(nNodes,nNodes);
			F_cols = sparse(nNodes,nNodes);
			F_diag = sparse(nNodes,nNodes);

			% copy replica rows of tensor to corresponding free rows
			F_rows(P.free.edge,:) = S(P.replica.edge,:);
			if length(P.free.corner) > 0
				F_rows(P.free.corner,:) = sum(S(P.replica.corner,:),1);
			end
			
			% copy replica cols of tensor to corresponding free cols
			F_cols(:,P.free.edge) = S(:,P.replica.edge);
			if length(P.free.corner) > 0
				F_cols(:,P.free.corner) = sum(S(:,P.replica.corner),2);
			end
			
			% copy replica diag entries of tensor to corresponding free diag entries
			F_diag(P.free.edge,:) = F_cols(P.replica.edge,:);
			if length(P.free.corner) > 0
				F_diag(P.free.corner,:) = sum(F_cols(P.replica.corner,:),1);
			end

			% assemble final result
			S_corr = F_rows + F_cols + F_diag;

		end

		function b_corr = computePeriodicCorrectionVector(self,b)

			% unpack variables
			nNodes = self.domain.mesh.nNodes;
			P = self.domain.boundary.P_nodes;

			% compute source correction vector
			b_corr = sparse(nNodes,1);
			b_corr(P.free.edge) = b(P.replica.edge);
			b_corr(P.free.corner) = sum(b(P.replica.corner));

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

			%{
			% copy solution to periodic replica nodes
			P = self.domain.boundary.P_nodes;
			self.solution(P.replica.edge,:) = self.solution(P.free.edge,:);
			self.solution(P.replica.corner,:) = ...
							repmat(self.solution(P.free.corner,:),3,1);					
			%self.solution(P.replica.corner,:) = self.solution(P.free.corner,:);
			%}

            % ensure solution is full
            self.solution = full(self.solution);

		end

		function [f,t,U,V] = checkVariables(self,f)

			% unpack variables
			nNodes    = self.domain.mesh.nNodes;

			% store function type data
			isCoupled     = num2str(Coefficients.isCoupled(f));
			isNonlinear   = num2str(Coefficients.isNonlinear(f));
			isTimeVarying = num2str(Coefficients.isTimeVarying(f));
			code = strcat(isCoupled,isNonlinear,isTimeVarying);
			code = str2num(code);

			% Update function signature and store variables and/or set dummy variables
			switch code
				% Spatially varying only
				case 000
					f = @(x1,x2,t,u,v)(f(x1,x2));
					t = 0;
					U = zeros(nNodes,1);
					V = 0;

				% Time varying
				case 001
					f = @(x1,x2,t,u,v)(f(x1,x2,t));
					t = self.t;
					U = zeros(nNodes,1);
					V = 0;
				
				% Nonlinear
				case 010
					f = @(x1,x2,t,u,v)(f(x1,x2,u));
					t = 0;
					U = self.U;
					V = 0;
				
				% Time varying and nonlinear
				case 011
					f = @(x1,x2,t,u,v)(f(x1,x2,t,u));
					t = self.t;
					U = self.U;
					V = 0;
				
				% Coupled
				case 100
					f = @(x1,x2,t,u,v)(f(x1,x2,v));
					t = 0;
					U = zeros(nNodes,1);
					%V = self.V;		% <~~~ use if V is not time-lagged
					V = self.V_prev;    % <~~~ use if V is time-lagged

				% Coupled and time varying
				case 101
					f = @(x1,x2,t,u,v)(f(x1,x2,t,v));
					t = self.t;
					U = zeros(nNodes,1);
					%V = self.V;		% <~~~ use if V is not time-lagged
					V = self.V_prev;    % <~~~ use if V is time-lagged

				% Coupled and nonlinear
				case 110
					f = @(x1,x2,t,u,v)(f(x1,x2,u,v));
					t = 0;
					U = self.U;
					%V = self.V;		% <~~~ use if V is not time-lagged
					V = self.V_prev;    % <~~~ use if V is time-lagged

				% Coupled, time varying, and nonlinear
				case 111
					f = @(x1,x2,t,u,v)(f(x1,x2,t,u,v));
					t = self.t;
					U = self.U;
					%V = self.V;		% <~~~ use if V is not time-lagged
					V = self.V_prev;    % <~~~ use if V is time-lagged
			end
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
			h = trisurf(elements3,coordinates(:,1),coordinates(:,2),U', ...
				'facecolor','interp');

			% format plot
			view(10,40);
			title('$u$','Interpreter','latex');

			% set zlim
			zMin = min(self.solution(:));
			zMax = max(self.solution(:));
			if zMin == zMax
				zMax = zMin + 1;
			end
			zlim([zMin, zMax]);

			% set xlim
			xMin = min(self.domain.mesh.nodes(:,1));
			xMax = max(self.domain.mesh.nodes(:,1));
			xlim([xMin, xMax]);

			% set ylim
			yMin = min(self.domain.mesh.nodes(:,2));
			yMax = max(self.domain.mesh.nodes(:,2));
			ylim([yMin, yMax]);

			% set colorbar
			colorbar;
			clim([zMin, zMax]);

			% adjust position 
			f = gcf;
			heightRatio = (yMax - yMin) / (xMax - xMin);
			width = 500;
			f.Position = [100, 100, width, width * heightRatio];

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
			segIDs2
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




		%{
		% UNVECTORIZED ASSEMBLY FUNCTIONS (FOR REFERENCE) --------------------%
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

		function area = elementArea(self,vertices)
			
			d = size(vertices,2);
			area = det([ones(1,d+1);vertices']) / prod(1:d);

		end

		% --------------------------------------------------------------------%
		%}


		%{
		% VECTORIZED ASSEMBLY FUNCTIONS (Attempt 1) --------------------------%
		% NOTE: In this model, the local tensors are computed as 3 x 3 x nElem
		% arrays, where each slice corresponds to the local tensor for a single
		% element. These local tensor arrays are summed to get an overall set of
		% local tensors called S_vec. This is scaled by the areas of the
		% elements, which completes centroid quadrature for the local tensors.
		% Then these local tensors are summed in a vectorized way to get the
		% global tensor S. This idea is technically faster than the method I
		% settled on above. But the run-time savings is very small; and at the
		% same time this approach requires fundamentally rewriting the way that
		% each solver comptues its solution. That is because currently, the
		% solvers compute global tensors individually, then these global tensors
		% are summed; but the approach below would require computing local
		% tensors, then summing, then making one global tensor. Ultimately, this
		% would require extensive rewrites to core parts of the code. At the
		% same time, testing shows that the approach below is only a few
		% hundredths of a second faster. So it's just not worth it. That said,
		% maybe in the future I will revisit this approach, so I am keeping it
		% here for reference.
		
		function S_vec = computeLocalTensors_vec(self)
		% NOTE: To implement this approach, each individual solver would have
		% its own version of the "computeLocalTensors_vec" function. Currently,
		% each solver has its own "pattern" for how it sums the global tensors.
		% That same pattern would appear in the computation of the local
		% tensors, obtaining some 3 x 3 x nElem output called S_vec, which would
		% then be assembled into the single global tensor S. Based on this
		% model, not only is the code below specific to the Poisson problem, but
		% the code below is incomplete since it does not account for the tensor
		% that is introduced by the Robin boundary values. You could (1) add
		% that tensor to the code below, or (2) add a separate function that
		% computes the Robin tensor, and then sum that tensor with the tensor
		% below. Either way, you can see that this approach requires extensive
		% rewrites to the code base.
		
			% store variables
			r = self.coefficients.r;

			% compute local tensors
			A_vec = self.localStiffnessMatrix_vec;
			M_vec = self.localMassMatrix_vec(r);

			% sum local matrices
			S_vec = A_vec + M_vec;

		end


		function self = assembleTensors_vec(self)
		% This function creates the local stiffness matrices S_vec, and then
		% assembles them into the global stiffness matrix S. The code below does
		% get used in the final approach above. It's just that I'm using
		% something like the code below to assemble each individual tensor into
		% its respective global counterpart. For example, A_vec becomes A; M_vec
		% becomes M, etc. Then I am summing these individual global tensors to
		% get the final global tensor S. 
		
			% get local tensors
			% note: computeLocalTensors_vec function is specific to each individual problem
			S_vec = self.computeLocalTensors_vec;

			% unpack variables
			nNodes     = self.domain.mesh.nNodes;
			elementInd = self.domain.mesh.elements(self.domain.mesh.effectiveElems,:);

			% scale S_vec by areas to complete centroid quadrature
			areas = self.domain.mesh.areas(self.domain.mesh.effectiveElems);
			areas = permute(areas,[3,2,1]);
			S_vec = pagemtimes(areas,S_vec);

			% Get summation locations
			[I, J] = ndgrid(1:3, 1:3);
			ind_vec = sub2ind([nNodes,nNodes], elementInd(:,I(:)), elementInd(:,J(:)));
			ind_vec = reshape(ind_vec',[],1);

			% sum the values
			% note: flag 'true' below indicates sparse matrix is to be created
			S_vec_flat = S_vec(:);
			S = accumarray(ind_vec,S_vec_flat,[nNodes * nNodes,1],[],[],true);
			S = reshape(S,[nNodes,nNodes]);

			% temporary: store matrices
			self.tensors.A = S;
			temp = size(self.tensors.A);
			self.tensors.M_r = sparse(temp(1),temp(2));			

		end
		%}

