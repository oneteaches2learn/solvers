classdef Domain2d

	properties
		domainArea
		mesh
		boundary
		time
	end

	methods
		% CONSTRUCTOR
		function self = Domain2d(x,y)

			if nargin == 2
				% set deconstructed geometry description matrix
				gd = Domain2d.getGeometryDescriptionMatrix(x,y);
				ns = Domain2d.getNameSpace;
				sf = Domain2d.getSetFunction;
				dl_mat = decsg(gd,sf,ns);

				% create boundary object
				self.boundary = Boundary2d(dl_mat);
			end

		end

		% SETTERS
		function self = setTime(self,T,dt,eq)

			if nargin == 2
				self.time = TimeStepping(T);
			elseif nagin == 3
				self.time = TimeStepping(T,dt);
			elseif nargin == 4
				self.time = TimeStepping(T,dt,eq);
			end

		end

		function self = setBoundary(self,bcTypes,bcConds,mesh)

			self = self.setBCTypes(bcTypes);
			self = self.setBCConditions(bcConds);
			
			if nargin == 4
				self = self.setBoundaryNodes(mesh);
			elseif nargin == 3 && ~isempty(self.mesh)
				self = self.setBoundaryNodes;
			end
		end

		function self = setBCTypes(self,bcTypes)
			
			self.boundary = self.boundary.setBCTypes(bcTypes);

		end

		function self = setBCConditions(self,bcConds)

			self.boundary = self.boundary.setBCConds(bcConds);

		end

		function self = setBoundaryNodes(self,mesh)

			if nargin == 1, mesh = self.mesh; end

			self.boundary = self.boundary.setEdgeNodes(mesh);
			self.boundary = self.boundary.setBoundaryNodeLists(mesh);
			self.boundary = self.boundary.setFreeNodes(mesh);

		end

		function self = setMesh(self,p,base,region)

			% region is unused for unpunctured domains
			if nargin == 2, region = 'skip'; end

			% generate the mesh
			self.mesh = Mesh2d(self.boundary.dl.mat,p,base);

			% set effective nodes and elements
			self = self.setEffectiveNodes;
			self = self.setEffectiveElements;

			% compute domain area
			self.domainArea = sum(self.mesh.areas);

		end

		function self = setEffectiveNodes(self)

			% store node data
			self.mesh.effectiveNodes = [1:1:self.mesh.nNodes];
			self.mesh.unusedNodes = [];

		end

		function self = setEffectiveElements(self)

			% store element data
			self.mesh.effectiveElems = [1:1:self.mesh.nElems];
			self.mesh.unusedElems = [];

		end


		% GETTERS
		function output = getElementCoordinates(self,nElem,requestedCoord)

			% store domain information
			coordinates = self.mesh.nodes;
			elements3 = self.mesh.elements;

			% store coordinates
			xElem = coordinates(elements3(nElem,:),1);
			yElem = coordinates(elements3(nElem,:),2);

			% set output
			if nargin == 3
				if requestedCoord == 1
					output = xElem;
				elseif requestedCoord == 2
					output = yElem;
				end
			else
				output = [xElem,yElem];
			end

		end

		function A = getDomainArea(self)

			% get domain area
			A = area(self.mesh.Mesh);

		end

		function unusedNodes = get_unusedNodes(self)

			nodeList = [1:1:self.mesh.nNodes];
			unusedNodes = setdiff(nodeList,self.mesh.effectiveNodes);

		end

		function unusedElems = get_unusedElems(self)

			elemList = [1:1:self.mesh.nElems];
			unusedElems = setdiff(elemList,self.mesh.effectiveElems);

		end


		% Y-LINE
		function self = add_yline(self,varargin)

			self.boundary = self.boundary.add_yline(varargin{:});
			if ~isempty(self.mesh)
				self = self.setMesh(self.mesh.p,self.mesh.base);
			end

		end


		% INCLUSIONS ON/OFF
		function self = inclusionsON(self,varargin)

			self.boundary = self.boundary.inclusionsON;
			if ~isempty(self.mesh)
				self = self.setMesh(self.mesh.p,self.mesh.base);
			end

		end

		function self = inclusionsOFF(self,varargin)

			self.boundary = self.boundary.inclusionsON;
			if ~isempty(self.mesh)
				self = self.setMesh(self.mesh.p,self.mesh.base);
			end

		end



		% UTILITY FUNCTIONS
		function vals = function2midpoints(self,varargin)
		% function2midpoints evaluates arguments on element midpoints
		% Inputs may be symfun, sym, function_handle, or vector. If sym, symfun,
		% or function_handle, the function value on the midpoints of the mesh
		% will be computed and stored as a vector. If input is vector, it should
		% be a vector of values computed on the mesh edge midpoints. Note that
		% vector arguments will be passed unchanged. 
		% 
		% For stationary problems, output will be a matrix of column vectors
		% where vals(:,i) corresponds to input varargin{i} on the mesh
		% midpoints. For time-varying problems, output will be a 3d array where
		% vals(:,i,:) corresponds to varagin{i} on the mesh midpoints with the
		% 3rd dimension corresponding to time.

			% store copy of inputs
			temp = varargin;

			% convert symfun to function_handle
			for i = 1:length(temp)
				if isa(temp{i},'symfun') || isa(temp{i},'sym')
					x = sym('x',[1 2]);
					if ~isempty(self.time)
						syms t;
						f = matlabFunction(symfun(temp{i},[x t]));
						temp{i} = @(x1,x2,t)(f(x1,x2,t) + zeros(length(x1),length(t)));
						% NOTE: the above function has been "wrapped" to ensure
						% that the output is a matrix with the correct
						% dimensions. This is necessary in the case that the
						% function is constant in space or in time (or both), in
						% which case MATLAB will output a vector (or scalar)
						% instead of a matrix. The wrapper is the `+
						% zeros(length(x1),length(t))` part of the function
						% definition, and this ensures that the output is a
						% matrix of the appropriate size.
						else
						f = matlabFunction(symfun(temp{i},x));
						temp{i} = @(x1,x2)(f(x1,x2) + zeros(length(x1),1));
					end
				end
			end

			% compute midpoint values
			for i = 1:length(temp)

				% if function_handle, compute values directly
				if isa(temp{i},'function_handle')
					func = temp{i};

					% for time-varying problems
					if nargin(temp{i}) == 3

						tGrid = self.time.tGrid;
						temp{i} = func(self.mesh.midpoints(:,1), ...
											self.mesh.midpoints(:,2),tGrid);

						% NOTE: The following two blocks of code were written
						% before the use of the wrapper function. I believe they
						% are no longer necessary, but are kept here in case
						% future issues arise with the wrapper function.
						%{
						% if uTrue is constant in space, then must be manually duplicated
						if size(temp{i},1) == 1
							sol = temp{i};
							for i = 1:size(sol,2)
								sol2(:,i) = sol(i) * ones(self.mesh.nEdges,1);
							end
							temp{i} = sol2;
						end
						%}

						%{
						% if uTrue is constant in time, then N_t-many copies must be manually made
						if size(temp{i},2) == 1
							temp{i} = repmat(temp{i},1,self.time.N_t);
						end
						%}

					% for stationary problems
					elseif nargin(temp{i}) == 2
						temp{i} = func(self.mesh.midpoints(:,1), ...
												self.mesh.midpoints(:,2));
					end

				% if vector, pass through
				elseif isa(temp{i},'double')
					...
				end

			end

			% store output as columns of non-sparse matrix
			for i = 1:length(temp)
				if ~isempty(self.time)
					vals(:,i,:) = temp{i};
				else
					vals(:,i) = temp{i};
				end
			end
			vals = full(vals);
		end
	
		function vals = function2centroids(self,varargin)
		% function2centroids evaluates arguments on element centroids. 
		% Inputs may be symfun, sym, function_handle, or vector. If sym, symfun,
		% or function_handle, the function value on the midpoints of the mesh
		% will be computed and stored as a vector. If input is vector, it should
		% be a vector of values computed on element centroids. Note that vector
		% arguments will be passed unchanged.
		% 
		% For stationary problems, output will be a matrix of column vectors
		% where vals(:,i) corresponds to input varargin{i} on the element 
		% centroids. For time-varying problems, output will be a 3d array where
		% vals(:,i,:) corresponds to varagin{i} on the element centroids with
		% the 3rd dimension corresponding to time.

			% store copy of inputs
			temp = varargin;

			% convert symfun to function_handle
			for i = 1:length(temp)
				if isa(temp{i},'symfun') || isa(temp{i},'sym')
					x = sym('x',[1 2]);
					if ~isempty(self.time)
						syms t;
						temp{i} = matlabFunction(symfun(temp{i},[x t]));
					else
						temp{i} = matlabFunction(symfun(temp{i},x));
					end
				end
			end

			% compute centroid values
			for i = 1:length(temp)

				% if function_handle, compute values directly
				if isa(temp{i},'function_handle')
					func = temp{i};

					% for time-varying problems
					if nargin(temp{i}) == 3

						tGrid = self.time.tGrid;
						temp{i} = func(self.mesh.centroids(:,1), ...
										self.mesh.centroids(:,2),tGrid);

						% if uTrue is constant in time, then N_t-many copies must be manually made
						if size(temp{i},2) == 1
							temp{i} = repmat(temp{i},1,self.time.N_t);
						end

					% for stationary problems
					elseif nargin(temp{i}) == 2
						temp{i} = func(self.mesh.centroids(:,1), ...
											self.mesh.centroids(:,2));
					end

				% if vector, pass through
				elseif isa(temp{i},'double')
					...
				end

			end

			% store output as columns of non-sparse matrix
			for i = 1:length(temp)
				if ~isempty(self.time)
					vals(:,i,:) = temp{i};
				else
					vals(:,i) = temp{i};
				end
			end
			vals = full(vals);

		end
	
		function vals = function2nodes(self,varargin)
		% function2nodes evaluates arguments on mesh nodes. 
		% Inputs may be symfun, sym, function_handle, or vector. If sym, symfun,
		% or function_handle, the function value on the mesh nodes will be
		% computed and stored as a vector. If input is vector, it should
		% be a vector of values computed on mesh nodes. Note that vector
		% arguments will be passed unchanged.
		% 
		% For stationary problems, output will be a matrix of column vectors
		% where vals(:,i) corresponds to input varargin{i} on the mesh
		% nodes. For time-varying problems, output will be a 3d array where
		% vals(:,i,:) corresponds to varagin{i} on the mesh nodes with the
		% 3rd dimension corresponding to time.

			% store copy of inputs
			temp = varargin;

			% convert symfun to function_handle
			for i = 1:length(temp)
				if isa(temp{i},'symfun') || isa(temp{i},'sym')
					x = sym('x',[1 2]);
					if ~isempty(self.time)
						syms t;
						temp{i} = matlabFunction(symfun(temp{i},[x t]));
					else
						temp{i} = matlabFunction(symfun(temp{i},x));
					end
				end
			end

			% compute nodal values
			for i = 1:length(temp)

				% if function_handle, compute values directly
				if isa(temp{i},'function_handle')
					func = temp{i};

					% for time-varying problems
					if nargin(temp{i}) == 3

						tGrid = self.time.tGrid;
						temp{i} = func(self.mesh.nodes(:,1), ...
										self.mesh.nodes(:,2),tGrid);

						% if uTrue is constant in time, then N_t-many copies must be manually made
						if size(temp{i},2) == 1
							temp{i} = repmat(temp{i},1,self.time.N_t);
						end

					% for stationary problems
					elseif nargin(temp{i}) == 2
						temp{i} = func(self.mesh.nodes(:,1), ...
											self.mesh.nodes(:,2));
					end

				% if vector, pass through
				elseif isa(temp{i},'double')
					...
				end

			end

			% store output as columns of non-sparse matrix
			for i = 1:length(temp)
				if ~isempty(self.time)
					vals(:,i,:) = temp{i};
				else
					vals(:,i) = temp{i};
				end
			end
			vals = full(vals);

		end

		function vals = nodes2midpoints(self,varargin)
		% nodes2midpoints converts nodal values to values on mesh midpoints
		% Inputs may be symfun, sym, function_handle, or vector. If sym, symfun,
		% or function_handle, the function value on the midpoints of the mesh
		% will be computed and stored as a vector. If input is vector, it should
		% be a vector of values computed on the mesh nodes. From this vector of
		% nodal values, function values on midpoints will be computed as
		% averages of nodal values. 
		%
		% For stationary problems, output will be a matrix of column vectors
		% where vals(:,i) corresponds to input varargin{i} on the mesh
		% midpoints. For time-varying problems, output will be a 3d array where
		% vals(:,i,:) corresponds to varagin{i} on the mesh midpoints with the
		% 3rd dimension corresponding to time.
		
			% store copy of inputs
			temp = varargin;

			% compute midpoint values
			for i = 1:length(temp)

				% if vector, take average of nodal values along edges
				if isa(temp{i},'double')
					vec = temp{i};
					temp{i} = (vec(self.mesh.edges(:,1),:) + vec(self.mesh.edges(:,2),:)) / 2;
				end

			end

			% compute symfun or function_handle values on midpoints
			vals = full(self.function2midpoints(temp{:}));

		end

		function vals = nodes2centroids(self,varargin)
		% nodes2centroids converts nodal values to values on element centroids
		% Inputs may be symfun, sym, function_handle, or vector. If sym, symfun,
		% or function_handle, the function value on the centroids of the mesh
		% will be computed and stored as a vector. If input is vector, it should
		% be a vector of values computed on the mesh nodes. From this vector of
		% nodal values, function values on centroids will be computed as
		% averages of nodal values. 
		%
		% For stationary problems, output will be a matrix of column vectors
		% where vals(:,i) corresponds to input varargin{i} on the mesh
		% midpoints. For time-varying problems, output will be a 3d array where
		% vals(:,i,:) corresponds to varagin{i} on the mesh midpoints with the
		% 3rd dimension corresponding to time.
		
			% store copy of inputs
			temp = varargin;

			% compute centroid values
			for i = 1:length(temp)

				% if vector, take average of nodal values for each element
				if isa(temp{i},'double')
					vec = temp{i};
					temp{i} = (vec(self.mesh.elements(:,1),:) ...
							 	+ vec(self.mesh.elements(:,2),:) ...
							 	+ vec(self.mesh.elements(:,3),:)) / 3;
				end

			end

			% compute symfun or function_handle values on midpoints
			vals = full(self.function2centroids(temp{:}));

		end


		function df = gradient_nodal(self,f)

			% the i-th page gives the coordinates of the i-th element
			X = self.mesh.elementCoords;
			X = rot90(X);
			X = permute(X,[1,3,2]);

			% take the difference of the coordinates
			A1 = X(2,:,:) - X(1,:,:);
			A2 = X(3,:,:) - X(1,:,:);
			A =  [A1; A2];

			% the i-th page gives the values of f on the i-th element
			Z = f(self.mesh.elements);
			Z = rot90(Z);
			Z = permute(Z,[1,3,2]);

			% take the difference of the values
			B1 = Z(2,:,:) - Z(1,:,:);
			B2 = Z(3,:,:) - Z(1,:,:);
			B = [B1; B2];

			% the i-th page gives the i-th gradient vector
			df = pagemldivide(A,B);

			% reshape, the i-th row gives the gradient of f on the i-th element
			df = permute(df,[3, 1, 2]);

			%{
			% below is a loop-based version of the code to illustrate how the
			% above vectorized code works, and to demonstrate that the
			% vectorized code gives the same results as the simlper loop-based
			% code
			for i = 1:self.mesh.nElems

				% get coordinates of element nodes
				Xbar(3,:) = self.mesh.nodes(self.mesh.elements(i,1),:);
				Xbar(2,:) = self.mesh.nodes(self.mesh.elements(i,2),:);
				Xbar(1,:) = self.mesh.nodes(self.mesh.elements(i,3),:);
				%Xbar - X(:,:,i)

				% get corresponding solution values
				zbar = f(self.mesh.elements(i,:));
				%zbar = flipud(zbar)

				% take difference of coordinates
				Abar = [Xbar(2,:) - Xbar(1,:); Xbar(3,:) - Xbar(1,:)];
				%A(:,:,i) - Abar

				% take difference of values
				Bbar = [zbar(2) - zbar(1); zbar(3) - zbar(1)];
				%B(:,:,i) - Bbar

				% compute gradient
				dfBar = Abar\Bbar;
				dfBar = dfBar';
				df = dfBar;
				%dfBar - df(i,:)

			end
			%}

		end

		function obj = refineMesh(self,numRefinements)
		% refineMesh(self) uses MATLAB's inbuilt tools to refine the mesh using
		% the 4-triangle refinement.

			if self.boundary.nEdges > 4
				error("Mesh refinement is not supported for domains with more than 4 edges.")
			end
		
			if nargin == 1
				numRefinements = 1;
			end

			% get inputs for refinement tool
			obj = self;
			[p,e,t] = meshToPet(obj.mesh.Mesh);
			dl_mat = obj.boundary.dl.mat;

			% sequentially refine the mesh
			for i = 1:numRefinements

				[p,e,t] = refinemesh(dl_mat,p,e,t);

			end

			% make dummy model and generate FEmesh
			model = createpde();
			geometryFromMesh(model,p,t(1:3,:));
			FEmsh = model.Mesh;

			% make Mesh2d object
			msh = Mesh2d(FEmsh);

			% update mesh
			obj.mesh = msh;

		end


		% PLOTTERS
		function h = plot(self,NameValueArgs)

			arguments
				self
				NameValueArgs.NodeLabels string = "off"
				NameValueArgs.NodeFontSize double = 22
				NameValueArgs.ElementLabels string = "off"
			end
			x = NameValueArgs;

			% Plot mesh
			hold on
			h = pdemesh(self.mesh.Mesh,ElementLabels=x.ElementLabels);

			% Label nodes, optional
			if x.NodeLabels == "on"
				nNodes = self.mesh.nNodes;	
				xData = self.mesh.nodes(:,1);
				yData = self.mesh.nodes(:,2);
				for i = 1:nNodes
					text(xData(i),yData(i),strcat('n',num2str(i)), ...
						'FontSize',x.NodeFontSize,'FontWeight','bold');
				end
			end

			hold off

		end

		function h = plotGeometry(self,NameValueArgs)

			arguments
				self
				NameValueArgs.EdgeLabels string = "off"
				NameValueArgs.FaceLabels string = "off"
			end
			x = NameValueArgs;

			h = pdegplot(self.boundary.dl.mat,FaceLabels=x.FaceLabels,EdgeLabels=x.EdgeLabels)

		end

		function x = x(self)
			x = [self.Vertices(1,1) self.Vertices(2,1)];
		end

		function y = y(self)
			y = [self.Vertices(1,2) self.Vertices(3,2)];
		end

		
		% THREE POINT QUADRATURE
		function int = threePointQuadrature(self,arg)
		% Computes Gaussian three point quadrature on elements of the mesh.
		% Input arg should be either a function_handle, a symbolic function, or
		% a vector of function evaluations on the mesh midpoints. If:
		%
		%	(1) arg is a function_handle, then arg will be evaluated on the
		%		midpoints of the edges of the mesh.
		%	(2) arg is a symbolic function, then ag will be converted to a
		%		function_handle and evaluated on the midpoints of the edges of
		%		the mesh.
		%	(3) arg is a vector, then arg should be indexed so that arg(i) stores the
		%		function evaluation on the i-th edge of the mesh, where the
		%		edges are ordered lexicographically. 
		%
		% Note: Gaussian three point quadrature has order four convergence as
		% the mesh is refined. It is accurate to machine precision for
		% polynomials up to order 2. 

			% process input into vector
			Arg = self.function2midpoints(arg);

			% compute quadrature
			elemEdges = self.mesh.elementEdges;
			elemAvg = sum([Arg(elemEdges(:,1),:,:), ...
							Arg(elemEdges(:,2),:,:),...
							Arg(elemEdges(:,3),:,:)],2) / 3;
			int = sum(self.mesh.areas .* elemAvg,"omitnan");
			int = reshape(int,[],1,1);

		end

		function IP = L2_IP_threePointQuadrature(self,arg1,arg2)
		% Computes the L2 inner product using Gaussian three point quadrature.
		% Inputs arg1,arg2 should be either function_handles, symbolic functions, or
		% vectors of function evaluations on the midpoints of the edges of the
		% mesh. A mixture of these is also allowed. See threePointQuadrature.

			% process inputs into vectors
			vals = self.function2midpoints(arg1,arg2);
			Arg1 = vals(:,1,:);
			Arg2 = vals(:,2,:);

			% Compute the inner product
			IP = self.threePointQuadrature(Arg1 .* Arg2);

		end	

		function norm = L2norm_threePointQuadrature(self,arg)
		% Computes the L2 norm using Gaussian three point quadrature.
		% Input arg should be either function_handle, symbolic function, or
		% vector of function evaluations on the midpoints of the edges of the
		% mesh. See threePointQuadrature.

			% process input into vector
			Arg = self.function2midpoints(arg);

			% compute the norm
			norm_squared = self.L2_IP_threePointQuadrature(Arg,Arg);
			norm = sqrt(norm_squared);

		end

		function err = L2err_threePointQuadrature(self,arg1,arg2)
		% Computes the L2 error using Gaussian three point quadrature.
		% Inputs arg1,arg2 should be either function_handles, symbolic functions, or
		% vectors of function evaluations on the midpoints of the edges of the
		% mesh. A mixture of these is also allowed. See threePointQuadrature.

			% process inputs into vectors
			vals = self.function2midpoints(arg1,arg2);
			Arg1 = vals(:,1,:);
			Arg2 = vals(:,2,:);

			% compute error
			err = self.L2norm_threePointQuadrature(Arg1 - Arg2);

		end

		function IP = L2_IP_threePointQuadrature_nodal(self,arg1,arg2)
		% Computes the L2 inner product using Gaussian three point quadrature.
		% Inputs arg1,arg2 should be either function_handles, symbolic functions, or
		% vectors of function evaluations on the nodes of the mesh. A mixture of
		% these is also allowed.
		%
		% NOTE: For domain Omega with triangulation T, this function is intended
		% to compute the L2 inner product of a piecewise linear function on T
		% and a function_handle or symfun on Omega. In particular, suppose arg1
		% is a vector of values on the nodes of T. Then arg1 is interpreted as a
		% piecewise linear function on T. If arg2 is a function_handle or
		% symfun, then arg2 is evaluated on the midpoints of Tgq. The L2 inner
		% product is then computed using Gaussian three point quadrature. 
		%
		% NOTE: If arg1 and arg2 are both vectors of nodal values, then this
		% function is equivalent to L2_IP_nodalQuadrature and the higher order
		% convergence is lost. 
		
			% convert nodal arguments to midpoint arguments
			vals = self.nodes2midpoints(arg1,arg2);
			Arg1 = vals(:,1,:);
			Arg2 = vals(:,2,:);

			% compute inner product
			IP = self.L2_IP_threePointQuadrature(Arg1,Arg2);

		end

		function err = L2err_threePointQuadrature_nodal(self,arg1,arg2)
		% Computes the L2 error using Gaussian three point quadrature.
		% Inputs arg1,arg2 should be either function_handles, symbolic functions, or
		% vectors of function evaluations on the nodes of the mesh. A mixture of
		% these is also allowed. See L2_IP_threePointQuadrature_nodal.

			% convert nodal arguments to midpoint arguments
			vals = self.nodes2midpoints(arg1,arg2);
			Arg1 = vals(:,1,:);
			Arg2 = vals(:,2,:);

			% compute inner product
			err = self.L2err_threePointQuadrature(Arg1,Arg2);

		end


		% CENTROID QUADRATURE
		function int = centroidQuadrature(self,arg)
		% Computes centroid quadrature on elements of the mesh.
		% Input arg should be either a function_handle, a symbolic function, or
		% a vector of function evaluations on the element centroids. If:
		%
		%	(1) arg is a function_handle, then arg will be evaluated on the
		%		midpoints of the element centroids.
		%	(2) arg is a symbolic function, then ag will be converted to a
		%		function_handle and evaluated on the element centroids.
		% 	(3) arg is a vector, then arg should be indexed so that arg(i) 
		% 		stores the function evaluation on the i-th element of the mesh.

			% process input into vector
			Arg = self.function2centroids(arg);

			% compute quadrature
			int = sum(self.mesh.areas .* Arg,"omitnan");
			int = reshape(int,[],1,1);

		end

		function IP = L2_IP_centroidQuadrature(self,arg1,arg2)
		% Computes the L2 inner product using centroid quadrature.
		% Inputs arg1,arg2 should be either function_handles, symbolic
		% functions, or vectors of function evaluations on element centroids. A
		% mixture of these is also allowed. See centroidQuadrature.
		
			% process inputs into vectors
			vals = self.function2centroids(arg1,arg2);
			Arg1 = vals(:,1,:);
			Arg2 = vals(:,2,:);

			% Compute the inner product
			IP = self.centroidQuadrature(Arg1 .* Arg2);

		end	

		function norm = L2norm_centroidQuadrature(self,arg)
		% Computes the L2 norm using centroid quadrature.
		% Input arg should be either function_handle, symbolic function, or
		% vector of function evaluations on the element centroids. See
		% centroidQuadrature.
		
			% process input into vector
			Arg = self.function2centroids(arg);

			% compute the norm
			norm_squared = self.L2_IP_centroidQuadrature(Arg,Arg);
			norm = sqrt(norm_squared);

		end

		function err = L2err_centroidQuadrature(self,arg1,arg2)
		% Computes the L2 error using centroid quadrature.
		% Inputs arg1,arg2 should be either function_handles, symbolic
		% functions, or vectors of function evaluations on the element
		% centroids. A mixture of these is also allowed. See centroidQuadrature.
		
			% process inputs into vectors
			vals = self.function2centroids(arg1,arg2);
			Arg1 = vals(:,1,:);
			Arg2 = vals(:,2,:);

			% compute error
			err = self.L2norm_centroidQuadrature(Arg1 - Arg2);

		end

		function IP = L2_IP_centroidQuadrature_nodal(self,arg1,arg2)
		% Computes the L2 inner product using centroid quadrature.
		% Inputs arg1,arg2 should be either function_handles, symbolic
		% functions, or vectors of function evaluations on the nodes of the
		% mesh. A mixture of these is also allowed.
		%
		% NOTE: For domain Omega with triangulation T, this function is intended
		% to compute the L2 inner product of a piecewise linear function on T
		% and a function_handle or symfun on Omega. In particular, suppose arg1
		% is a vector of values on the nodes of T. Then arg1 is interpreted as a
		% piecewise linear function on T. If arg2 is a function_handle or
		% symfun, then arg2 is evaluated on the centroids of T. The L2 inner
		% product is then computed using centroid quadrature. 
		%
		% NOTE: If arg1 and arg2 are both vectors of nodal values, then this
		% function is equivalent to L2_IP_nodalQuadrature and the higher order
		% convergence is lost. 
		
			% convert nodal arguments to midpoint arguments
			vals = self.nodes2centroids(arg1,arg2);
			Arg1 = vals(:,1,:);
			Arg2 = vals(:,2,:);

			% compute inner product
			IP = self.L2_IP_centroidQuadrature(Arg1,Arg2);

		end

		function err = L2err_centroidQuadrature_nodal(self,arg1,arg2)
		% Computes the L2 error using centroid quadrature.
		% Inputs arg1,arg2 should be either function_handles, symbolic
		% functions, or vectors of function evaluations on the nodes of the
		% mesh. A mixture of these is also allowed. See
		% L2_IP_centroidQuadrature_nodal.

			% convert nodal arguments to midpoint arguments
			vals = self.nodes2centroids(arg1,arg2);
			Arg1 = vals(:,1,:);
			Arg2 = vals(:,2,:);

			% compute inner product
			err = self.L2err_centroidQuadrature(Arg1,Arg2);

		end


		% NODAL QUADRATURE
		function [int,elemInt] = nodalQuadrature(self,arg)
		% Computes nodal quadrature on elements of the mesh.
		% Input arg should be either a function_handle, a symbolic function, or
		% a vector of function evaluations on the mesh nodes. If:
		%
		%	(1) arg is a function_handle, then arg will be evaluated on the
		%		mesh nodes.
		%	(2) arg is a symbolic function, then arg will be converted to a
		%		function_handle and evaluated on the mesh nodes.
		% 	(3) arg is a vector, then arg should be indexed so that arg(i) 
		% 		stores the function evaluation on the i-th node of the mesh.

			% process input into vector
			Arg = self.function2nodes(arg);

			% compute quadrature
			nodalAvg = sum([Arg(self.mesh.elements(:,1),:,:), ...
							Arg(self.mesh.elements(:,2),:,:), ...
							Arg(self.mesh.elements(:,3),:,:)],2) / 3;
			int = sum(self.mesh.areas .* nodalAvg,"omitnan");
			int = reshape(int,[],1,1);

			% store element integrals
			elemInt = nodalAvg .* self.mesh.areas;

		end

		function IP = L2_IP_nodalQuadrature(self,arg1,arg2)
		% Computes the L2 inner product using nodal quadrature.
		% Inputs arg1,arg2 should be either function_handles, symbolic
		% functions, or vectors of function evaluations on the mesh nodes
		% A mixture of these is also allowed. See nodalQuadrature.
		
			% process inputs into vectors
			vals = self.function2nodes(arg1,arg2);
			Arg1 = vals(:,1,:);
			Arg2 = vals(:,2,:);

			% Compute the inner product
			IP = self.nodalQuadrature(Arg1 .* Arg2);

		end	

		function norm = L2norm_nodalQuadrature(self,arg)
		% Computes the L2 norm using nodal quadrature.
		% Input arg should be either function_handle, symbolic function, or
		% vector of function evaluations on the mesh nodes. See nodalQuadrature.
		
			% process input into vector
			Arg = self.function2nodes(arg);

			% compute the norm
			norm_squared = self.L2_IP_nodalQuadrature(Arg,Arg);
			norm = sqrt(norm_squared);

		end

		function err = L2err_nodalQuadrature(self,arg1,arg2)
		% Computes the L2 error using nodal quadrature.
		% Inputs arg1,arg2 should be either function_handles, symbolic
		% functions, or vectors of function evaluations on the mesh nodes. See
		% nodalQuadrature.
		
			% process inputs into vectors
			vals = self.function2nodes(arg1,arg2);
			Arg1 = vals(:,1,:);
			Arg2 = vals(:,2,:);

			% compute error
			err = self.L2norm_nodalQuadrature(Arg1 - Arg2);

		end

	end

	methods (Static)
		
		function gd = getGeometryDescriptionMatrix(x,y)

			% geometry description matrix
			gd = zeros(10,1);
			gd(:,1) = [3 4 x(1) x(2) x(2) x(1) y(1) y(1) y(2) y(2)]';

		end

		function ns = getNameSpace()

			% name space
			ns = char('BD');
			ns = ns';

		end

		function sf = getSetFunction()

			% set function
			sf = 'BD';

		end

		function gd = gd_from_vertices(vertices)

			nSquares = size(vertices,1) / 4;
			gd = zeros(10,nSquares);

			for i = 0:nSquares-1
				
				% get x and y limits
				coord = vertices([1:4] + 4*i,:);
				xLim = [min(coord(:,1)) max(coord(:,1))];
				yLim = [min(coord(:,2)) max(coord(:,2))];

				% set gd column
				gd(:,i+1) = [3; 4; Domain2d.vertices_from_bounds(xLim,yLim)];

			end

		end

		function vert = vertices_from_bounds(xLim,yLim)

			vert = [xLim(1) xLim(2) xLim(2) xLim(1) yLim(1) yLim(1) yLim(2) yLim(2)]';

		end

		function dom = domainFromGmsh(msh)
		
			% STEP 2: store the elements; 
			%   trim off unnecessary fourth column (bc elements are 3d)
			%   transpose because matlab expects a 3 x nElem matrix
			%   reorder the nodes to match matlab's ordering
			elements = msh.TRIANGLES(:,1:3);
			elements = elements';
			%elements = elements([1 3 2],:);


			% STEP 3: store the nodes;
			%   trim off the third column (bc nodes are 2d)
			%   transpose because matlab expects a 2 x nNode matrix
			nodes = msh.POS(:,1:2);
			nodes = nodes';


			% STEP 4: create the mesh object
			%   to use matlab's pde toolbox, define a model
			%   the mesh and geometry will be stored in the model
			model = createpde();
			geo = geometryFromMesh(model,nodes,elements);
			%temp_mesh = generateMesh(model,Hmax=2^(-5),GeometricOrder='linear');


			% STEP 5: create Tyler's mesh object
			mesh = Mesh2d();
			mesh.base = 2;
			mesh.p = 4;
			mesh.h = 1 / mesh.base^mesh.p;
			%mesh.Mesh = temp_mesh;
			mesh.Mesh = model.Mesh;
			mesh.effectiveNodes = [1:1:mesh.nNodes];
			mesh.effectiveElems = [1:1:mesh.nElems];
			mesh.areas = mesh.get_areas();


			% STEP 6: create Tyler's boundary object
			bdry = Boundary2d();

			% get endpoints of edge 20 from model.geometry
			edges = zeros(model.Geometry.NumEdges,4);
			for i = 1:model.Geometry.NumEdges
				nodes_i = findNodes(model.Mesh,"region","Edge",i);
				vert_1 = model.Mesh.Nodes(:,nodes_i(1))';
				vert_2 = model.Mesh.Nodes(:,nodes_i(end))';
				edges(i,:) = [vert_1, vert_2];
			end
			[edges,edgeDict] = Domain2d.orderEdges(edges);

			% create edge list for boundary object
			%{
			% create one edge object for each edge
			edgesObj = [];
			for i = 1:size(edges,1)
				vert1 = edges(i,1:2);
				vert2 = edges(i,3:4);
				edge_i = BoundaryEdge2d(vert1,vert2);
				edge_i.segIDs = i;
				edge_i.nodes = findNodes(model.Mesh,"region","Edge",edgeDict(i));
				edgesObj = [edgesObj, edge_i];
			end
			%}

			% hack: create four total edges
			edgesObj = [];

			% first three edges are the first three edges
			for i = 1:3
				vert1 = edges(i,1:2);
				vert2 = edges(i,3:4);
				edge_i = BoundaryEdge2d(vert1,vert2);
				edge_i.segIDs = i;
				edge_i.nodes = findNodes(model.Mesh,"region","Edge",edgeDict(i));
				edgesObj = [edgesObj, edge_i];
			end

			% last edge contains all data from all remaining edges
			vert1 = edges(4,1:2);
			vert2 = edges(end,3:4);
			edge_i = BoundaryEdge2d(vert1,vert2);
			edge_i.segIDs = [4:1:size(edges,1)];
			edge_i.nodes = [];
			for i = 4:size(edges,1)
				edge_i.nodes = [edge_i.nodes, findNodes(model.Mesh,"region","Edge",edgeDict(i))];
			end
			edge_i.nodes = unique(edge_i.nodes);
			edgesObj = [edgesObj, edge_i];

			% add edges to boundary object
			bdry = Boundary2d();
			bdry.edges = edgesObj;
			bdry.gmshBoundary = 1;


			% STEP 7: create Tyler's domain object
			dom = Domain2d();
			dom.mesh = mesh;
			dom.boundary = bdry;
			dom.domainArea = sum(dom.mesh.areas);

		
		end


		function [ordered_edges,edgeDict] = orderEdges(edges)
		% Reorder edges to form a counterclockwise cycle around the boundary
		%
		% INPUTS:
		%   edges: n x 4 matrix, where each row is an edge defined by its two
		%          endpoints
		%
		% OUTPUTS:
		%   ordered_edges: m x 4 matrix, where each row is an edge defined by its
		%                  two endpoints, and the edges form a counterclockwise cycle
		%                  around the boundary
		%   edgeDict: m x 1 vector, where each element is the index of the edge in
		%             the original 'edges' matrix
			
			% STEP 1: find the starting vertex
			% extract vertex list from edge list
			vertices = unique([edges(:,1:2); edges(:,3:4)],'rows');
			
			% find the starting vertex (minimum y, then minimum x)
			min_y = min(vertices(:,2));
			min_y_points = vertices(vertices(:,2)==min_y,:);
			[min_x, idx_min_x] = min(min_y_points(:,1));
			startVertex = min_y_points(idx_min_x,:);
			

			% STEP 2: find the starting edge
			% find edges containing the starting vertex
			startEdges = find( (edges(:,1) == startVertex(1) & edges(:,2) == startVertex(2)) | ...
							(edges(:,3) == startVertex(1) & edges(:,4) == startVertex(2)) );
									
			% choose the edge representing counterclockwise traversal
			angles = zeros(length(startEdges),1);
			for i = 1:length(startEdges)

				% put startVertex at the begining of the edge
				idx = startEdges(i);
				edge = edges(idx,:);
				if edge(1:2) == startVertex
					otherVertex = edge(3:4);
				else
					otherVertex = edge(1:2);
					edges(idx,:) = [startVertex, otherVertex]; % Reverse edge
				end

				% calculate the angles these edges make with the x-axis
				delta = otherVertex - startVertex;
				angles(i) = atan2(delta(2), delta(1));
			end

			% choose the edge with the smallest angle to be the starting edge
			[~, minAngleIdx] = min(angles);
			startEdgeIdx = startEdges(minAngleIdx);
			startEdge = edges(startEdgeIdx,:);

			% update the ordered list and dictionary; mark the starting edge as used
			ordered_edges = startEdge;
			used_edges = false(size(edges,1),1);
			used_edges(startEdgeIdx) = true;
			edgeDict = [startEdgeIdx];
			
			
			% STEP 3: traverse the remaining edges
			currentVertex = startEdge(3:4);
			while sum(used_edges) < size(edges,1)

				% find the remaining edge that is connected to the current vertex
				connectedEdgeIdx = find(~used_edges & (...
					(edges(:,1) == currentVertex(1) & edges(:,2) == currentVertex(2)) | ...
					(edges(:,3) == currentVertex(1) & edges(:,4) == currentVertex(2))));
				if isempty(connectedEdgeIdx)
					break
				end

				% put vertices of connected edge in the correct order
				connectedEdge = edges(connectedEdgeIdx,:);
				if connectedEdge(1) == currentVertex(1) && connectedEdge(2) == currentVertex(2)
					nextVertex = connectedEdge(3:4);
				else
					nextVertex = connectedEdge(1:2);
					connectedEdge = [currentVertex, nextVertex]; % Reverse edge
				end

				% store the connected edge and mark it as used
				ordered_edges = [ordered_edges; connectedEdge];
				used_edges(connectedEdgeIdx) = true;
				edgeDict = [edgeDict; connectedEdgeIdx];

				% update the current vertex
				currentVertex = nextVertex;

			end

			% The 'ordered_edges' variable now contains the edges in counterclockwise order
		end
    

	end
end
