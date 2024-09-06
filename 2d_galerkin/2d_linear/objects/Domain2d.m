classdef Domain2d

	properties
		geometryMatrix
		dl
		base
		p
		h
		edges
		effectiveNodes
		effectiveElems
		unusedNodes
		unusedElems
		nEdges
		domainArea
		boundaryNodes
		freeNodes
		xLim
		yLim
		mesh
		time
	end

	methods
		% CONSTRUCTOR
		function self = Domain2d(x,y)

			% set deconstructed geometry description matrix
			gd = Domain2d.getGeometryDescriptionMatrix(x,y);
			ns = Domain2d.getNameSpace;
			sf = Domain2d.getSetFunction;
			self.dl = decsg(gd,sf,ns);

			% store geometry matrix
			self.geometryMatrix = gd;

			% store xLim and yLim
			self.xLim = x;
			self.yLim = y;

			% set and store edges
			self.edges = self.setEdgeGeometry;
			self.nEdges = self.get_nEdges;

		end

		% SETTERS
		function self = set.time(self,time)

			self.time = time;

		end

		function edges = setEdgeGeometry(self)

			% store variables
			dl = self.dl

			% detect columns of dl where x-coordinate is constant
			xConst = ((dl(2,:) - dl(3,:)) == 0);

			% mark which columns of dl designated the start of a new edge
			newEdge = [1, (xConst(1:end-1) - xConst(2:end)) ~= 0];

			% indicate which dl cols correspond to each edge 
			edgeNums(1) = 1;
			for i = 2:length(newEdge)
				edgeNums(i) = edgeNums(i-1) + newEdge(i);
			end

			% create edges
			edges = [];
			for i = 1:4

				% make map between geometry edges and dl columns
				IDs = find(edgeNums - i == 0);

				% get start and stop coordinates
				startCoord = [dl(2,IDs(1)),dl(4,IDs(1))];
				stopCoord = [dl(3,IDs(end)),dl(5,IDs(end))];

				% generate edges
				edge_i = BoundaryEdge2d(startCoord,stopCoord);
				edge_i.ID = IDs;
				edges = [edges edge_i];
			end
		end

		function edges = setEdgeGeometry_old(self)

			% store variables
			gd = self.geometryMatrix;

			% get vertices
			vert = zeros(4,2);
			for i = 1:4
				vert(i,1) = gd(i+2,1);
				vert(i,2) = gd(i+6,1);
			end
			vert(5,:) = vert(1,:);

			% get midpoints
			for i = 1:4
				mdpt(i,:) = [(vert(i,1) + vert(i+1,1))/2, (vert(i,2) + vert(i+1,2))/2];
			end

			% set edges
			edges = [];
			for i = 1:4
				edge_i = BoundaryEdge2d(vert(i,:),vert(i+1,:));
				edge_i.ID = i;
				edges = [edges edge_i];
			end

		end

		function self = setEdgeBCTypes(self,boundary)

			edges = self.edges;

			for i = 1:self.nEdges
				edges(i).boundaryType = boundary.boundaryTypes{i};
			end

			self.edges = edges;

		end

		function self = setEdgeBCConditions(self,boundary)

			edges = self.edges;

			for i = 1:self.nEdges
				edges(i).boundaryCondition = boundary.boundaryConditions{i};
			end

			self.edges = edges;

		end

		function self = setMesh(self,p,base)

			% store inputs
			self.p = p;
			self.base = base;

			% compute h
			self.h = base^-p;
			
			% generate the mesh
			self.mesh = self.generateMesh;

			% store node data
			self.effectiveNodes = [1:1:self.mesh.nNodes];
			self.unusedNodes = [];

			% store element data
			self.effectiveElems = [1:1:self.mesh.nElems];
			self.unusedElems = [];

			% distribute boundary nodes to edges
			self.edges = self.distributeBoundaryNodes;

			% compile boundary node lists
			self.boundaryNodes.D = [];
			self.boundaryNodes.N = [];
			self.boundaryNodes.R = [];
			for i = 1:length(self.edges)
				if strcmp(self.edges(i).boundaryType,'D')
					self.boundaryNodes.D = [self.boundaryNodes.D self.edges(i).nodes];
				elseif strcmp(self.edges(i).boundaryType,'N')
					self.boundaryNodes.N = [self.boundaryNodes.N self.edges(i).nodes];
				elseif strcmp(self.edges(i).boundaryType,'R')
					self.boundaryNodes.R = [self.boundaryNodes.R self.edges(i).nodes];
				end
			end

			% sort D nodes and remove duplicates
			self.boundaryNodes.D = unique(self.boundaryNodes.D);

			% remove D nodes from R and N node lists, sort and remove duplicates
			self.boundaryNodes.N = setdiff(self.boundaryNodes.N,self.boundaryNodes.D);
			self.boundaryNodes.R = setdiff(self.boundaryNodes.R,self.boundaryNodes.D);

			% store free nodes
			self.freeNodes = setdiff(1:self.mesh.nNodes,self.boundaryNodes.D);

		end

		% GETTERS
		function [bNodes, nearestEdge] = getBoundaryNodes(self)

			% instantiate storage
			bNodes = [];
			nearestEdge = [];

			% collect node pairs for each edge
			for edgeID = 1:self.nEdges
				edgeNodes = self.mesh.Mesh.findNodes('region','Edge',edgeID);
				for j = 1:length(edgeNodes)-1
					bNodes = [bNodes; edgeNodes(j) edgeNodes(j+1)];
					nearestEdge = [nearestEdge; edgeID];
				end
			end
		end

		function nEdges = get_nEdges(self)

			nEdges = length(self.edges);

		end

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
			unusedNodes = setdiff(nodeList,self.effectiveNodes);

		end

		function unusedElems = get_unusedElems(self)

			elemList = [1:1:self.mesh.nElems];
			unusedElems = setdiff(elemList,self.effectiveElems);

		end

		% UTILITY FUNCTIONS
		function edges = distributeBoundaryNodes(self)

			edges = self.edges;
			for i = 1:length(self.edges)
				edgeID = edges(i).ID;
				edges(i).nodes = self.mesh.Mesh.findNodes('region','Edge',edgeID);
				edges(i).nNodes = length(edges(i).nodes);
			end

		end

		function mesh = generateMesh(self)

			geo = fegeometry(self.dl);
			geo = geo.generateMesh(Hmax=self.h,GeometricOrder='linear');
			mesh = Mesh2d(geo.Mesh);

		end

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
						temp{i} = matlabFunction(symfun(temp{i},[x t]));
					else
						temp{i} = matlabFunction(symfun(temp{i},x));
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

						% if uTrue is constant in time, then N_t-many copies must be manually made
						if size(temp{i},2) == 1
							temp{i} = repmat(temp{i},1,self.time.N_t);
						end

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

		function self = add_domain_y_line(self,varargin)

			% store variables
			dl = self.dl;
			yBar = varargin{1};
			
			% split region
			for i = [size(dl,2):-1:1]

				% if y_bar is between y_start and y_stop, then split the region
				if (yBar > dl(4,i) && yBar < dl(5,i)) || ...
				   (yBar > dl(5,i) && yBar < dl(4,i))

					% duplicate the i-th column
					dl = [dl(:,1:i), dl(:,i:end)];
					
					% adjust y_start and y_stop
					dl(5,i) = yBar;
					dl(4,i+1) = yBar;
				end
			end

			% increment region ids for those regions above yBar
			increment_cols = find(sum(dl(4:5,:) > yBar,1) > 0);
			dl(6,increment_cols) = dl(6,increment_cols) + 1;
			
			% add dividing line
			dl = [dl,[2,self.xLim(1),self.xLim(2),yBar,yBar,2,1]'];

			% store result
			self.dl = dl;

			self.edges = self.setEdgeGeometry;

		end

		% PLOTTING FUNCTIONS
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

			h = pdegplot(self.dl,FaceLabels=x.FaceLabels,EdgeLabels=x.EdgeLabels)

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
		function int = nodalQuadrature(self,arg)
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

	end
end
