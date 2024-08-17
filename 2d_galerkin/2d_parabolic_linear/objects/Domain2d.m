classdef Domain2d < fegeometry

	properties
		geometryMatrix
		base
		p
		h
		edges
		nNodes
		nElem
		elemAreas
		elemEdges
		domainArea
		boundaryNodes
		freeNodes
	end

	methods
		% CONSTRUCTOR
		function self = Domain2d(x,y)

			gd = Domain2d.getGeometryDescriptionMatrix(x,y);
			ns = Domain2d.getNameSpace;
			sf = Domain2d.getSetFunction;

			% call fegeometry superclass constructor
			g = decsg(gd,sf,ns);
			self@fegeometry(g);

			% store geometry matrix
			self.geometryMatrix = gd;

			% set and store edges
			self.edges = self.setEdgeGeometry;

		end

		% SETTERS
		function edges = setEdgeGeometry(self)

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

			% get edge IDs
			for i = 1:4
				edge_id(i) = self.nearestEdge(mdpt(i,:));
			end

			% get normal vectors
			n = [0 -1; 1 0; 0 1; -1 0];

			% set edges
			edges = [];
			for i = 1:self.NumEdges
				edge_i = BoundaryEdge2d(vert(i,:),vert(i+1,:),n(i,:));
				edge_i.ID = edge_id(i);
				edges = [edges edge_i];
			end

		end

		function self = setEdgeBCTypes(self,boundary)

			edges = self.edges;

			for i = 1:self.NumEdges
				edges(i).boundaryType = boundary.boundaryTypes{i};
			end

			self.edges = edges;

		end

		function self = setEdgeBCConditions(self,boundary)

			edges = self.edges;

			for i = 1:self.NumEdges
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
			
			% generate mesh
			self = self.generateMesh(Hmax=self.h,GeometricOrder='linear');

			% store number of nodes
			self.nNodes = size(self.Mesh.Nodes,2);

			% store number of elements
			self.nElem = size(self.Mesh.Elements,2);

			% store domain and element areas
			[self.domainArea,self.elemAreas] = area(self.Mesh);

			% store element edge numbers
			self.elemEdges = self.getElementEdges;

			% distribute boundary nodes to edges
			for i = 1:self.NumEdges
				edgeID = self.edges(i).ID;
				self.edges(i).nodes = self.Mesh.findNodes('region','Edge',edgeID);
				self.edges(i).nNodes = length(self.edges(i).nodes);
			end

			% compile boundary node lists
			self.boundaryNodes.D = [];
			self.boundaryNodes.N = [];
			self.boundaryNodes.R = [];
			for i = 1:self.NumEdges
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
			self.freeNodes = setdiff(1:self.nNodes,self.boundaryNodes.D);

		end

		% GETTERS
		function [bNodes, nearestEdge] = getBoundaryNodes(self)

			% instantiate storage
			bNodes = [];
			nearestEdge = [];

			% collect node pairs for each edge
			for edgeID = 1:self.NumEdges
				edgeNodes = self.Mesh.findNodes('region','Edge',edgeID);
				for j = 1:length(edgeNodes)-1
					bNodes = [bNodes; edgeNodes(j) edgeNodes(j+1)];
					nearestEdge = [nearestEdge; edgeID];
				end
			end
		end

		function output = getElementCoordinates(self,nElem,requestedCoord)

			% store domain information
			coordinates = self.Mesh.Nodes';
			elements3 = self.Mesh.Elements';

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

		function output = getElementArea(self,nElem)

			% store domain information
			coordinates = self.Mesh.Nodes';
			elements    = self.Mesh.Elements';

			% store element information
			elementInd    = elements(nElem,:);
			elementCoord  = coordinates(elementInd,:);

			% compute area
			output = det([1,1,1; elementCoord']) / 2;

		end

		function A = getDomainArea(self)

			% get domain area
			A = area(self.Mesh);

		end

		function G = getMeshGraph(self)

			% get mesh edges
			edges = self.getMeshEdges;

			% generate graph
			G = graph(edges(:,1),edges(:,2));

		end			


		function edges = getMeshEdges(self)
		% get a list of all edges, ordered lexicographically


			% get list of all element edges
			elemNodesOrig = self.Mesh.Elements;
			elemNodesPerm = elemNodesOrig([2 3 1],:);
			r = reshape(elemNodesOrig,[],1);
			s = reshape(elemNodesPerm,[],1);
			edges = [r,s];

			% edge nodes should be listed in ascending order
			edges = sort(edges,2);

			% remove duplicate edges (also sorts lexicographically)
			edges = unique(edges,'rows');

			% store as 'large enough' integer type
			if self.nNodes > intmax("uint32")
				edges = uint64(edges);
			elseif self.nNodes > intmax("uint16")
				edges = uint32(edges);
			elseif self.nNodes > intmax("uint8")
				edges = uint16(edges);
			else
				edges = uint8(edges);
			end
	
		end

		function elemEdges = getElementEdges(self)
		% creates a matrix of edges where each row represents an element and
		% each entry in that row is the lexicographic number of an edge of that
		% element

			% get edge and element data
			edges     = self.getMeshEdges;
			elemNodes = sort(self.Mesh.Elements',2);
			
			% store individual edges by element
			edge1 = [elemNodes(:,1),elemNodes(:,2)];
			edge2 = [elemNodes(:,1),elemNodes(:,3)];
			edge3 = [elemNodes(:,2),elemNodes(:,3)];

			% store edge locations
			[lia,elemEdges(:,1)] = ismember(edge1,edges,"rows");
			[lia,elemEdges(:,2)] = ismember(edge2,edges,"rows");
			[lia,elemEdges(:,3)] = ismember(edge3,edges,"rows");

			% store as 'large enough' integer type
			if size(edges,1) > intmax("uint32")
				elemEdges = uint64(elemEdges);
			elseif size(edges,1) > intmax("uint16")
				elemEdges = uint32(elemEdges);
			elseif size(edges,1) > intmax("uint8")
				elemEdges = uint16(elemEdges);
			else
				elemEdges = uint8(elemEdges);
			end
		
		end 

		function midptCoords = getMeshEdgeMidpoints(self,edges)

			if nargin == 1
				edges = self.getMeshEdges;
			end

			% store variables
			elemCoord = self.Mesh.Nodes';

			% compute midpoints
			startCoords = elemCoord(edges(:,1),:);
			endCoords   = elemCoord(edges(:,2),:);
			midptCoords = (startCoords + endCoords) / 2;

		end

		function centroids = getCentroids(self)

				% get centroid coordinates
				elemNodes = self.Mesh.Elements';
				elemCoord = self.Mesh.Nodes';

				% isolate coordanates of each node for each element
				nodes1 = elemCoord(elemNodes(:,1),:);
				nodes2 = elemCoord(elemNodes(:,2),:);
				nodes3 = elemCoord(elemNodes(:,3),:);

				% compute average of element nodes for each element
				centroids = (nodes1 + nodes2 + nodes3) / 3;
				
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
			h = self.pdemesh( ...
					ElementLabels=x.ElementLabels);

			% Label nodes, optional
			if x.NodeLabels == "on"
				nNodes = size(self.Mesh.Nodes,2);	
				xData = self.Mesh.Nodes(1,:);
				yData = self.Mesh.Nodes(2,:);
				for i = 1:nNodes
					text(xData(i),yData(i),strcat('n',num2str(i)), ...
						'FontSize',x.NodeFontSize,'FontWeight','bold');
				end
			end

			hold off

			pause();
			close();
		end

		function h = plotGeometry(self,NameValueArgs)

			arguments
				self
				NameValueArgs.EdgeLabels string = "off"
				NameValueArgs.FaceLabels string = "off"
			end
			x = NameValueArgs;

			h = self.pdegplot(FaceLabels=x.FaceLabels,EdgeLabels=x.EdgeLabels)

		end

		function x = x(self)
			x = [self.Vertices(1,1) self.Vertices(2,1)];
		end

		function y = y(self)
			y = [self.Vertices(1,2) self.Vertices(3,2)];
		end

		
		% ANALYSIS FUNCTIONS
		function int = nodalQuadrature(self,f)
		% Computes nodal quadrature on elements of the mesh.
		% Input f should be either a function_handle, a symbolic function, or a
		% vector of function evaluations on the nodes of the mesh. If:
		%
		%	(1) f is a function_handle, then f will be evaluated on the
		%		nodes of the mesh.
		%	(2) f is a symbolic function, then f will be converted to a
		%		function_handle and evaluated on the nodes of the mesh.
		%	(3) f is a vector, then f should be indexed so that f(i) stores the
		%		function evaluation on the i-th node of the mesh, where the
		%		node ordering is determined by self.Mesh.Nodes
		%
		% Note: Nodal quadrature has order two convergence as the mesh is
		% refined. It is accurate to machine precision for order 1 polynomials. 
			
			% if f is a function
			if isa(f,'function_handle') || isa(f,'symfun')

				% convert symfun to function_handle
				if isa(f,'symfun'), f = matlabFunction(f); end

				% compute F on edge midpoints 
				F = f(self.Mesh.Nodes(1,:)',self.Mesh.Nodes(2,:)');

			% else if f is a vector of nodal values
			elseif isvector(f)

				F = f; 

				% convert to column, if necessary
				if isrow(F), F = F'; end 

			end

			% compute average nodal values on each element
			elemNodes = self.Mesh.Elements';
			elemAvg = sum([F(elemNodes(:,1)), ...
					F(elemNodes(:,2)),F(elemNodes(:,3))],2) / 3;

			% dot product with element areas
			int = sum(self.elemAreas' .* elemAvg);

		end

		function int = centroidQuadrature(self,f)
		% Computes centroid quadrature on elements of the mesh.
		% Input f should be either a function_handle, a symbolic function, or a
		% vector of function evaluations on the centroids of the mesh. If:
		%
		%	(1) f is a function_handle, then f will be evaluated on the
		%		centroids of the mesh.
		%	(2) f is a symbolic function, then f will be converted to a
		%		function_handle and evaluated on the centroids of the mesh.
		%	(3) f is a vector, then f should be indexed so that f(i) stores the
		%		function evaluation on the centroid of the i-th element of the
		%		mesh, where the ordering of the elements is the same as in
		%		self.Mesh.Elements
		%
		% Note: Centroid quadrature has order two convergence as the mesh is
		% refined. It is accurate to machine precision for order 1 polynomials. 
			
			if isa(f,'function_handle') || isa(f,'symfun')

				% convert symfun to function_handle
				if isa(f,'symfun'), f = matlabFunction(f); end
				
				% compute F on centroids
				centroids = self.getCentroids;
				F = f(centroids(:,1),centroids(:,2));

			% else if f is a vector of centroid values
			elseif isvector(f)

				F = f; 

				% convert to column, if necessary
				if isrow(F), F = F'; end 

			end
			% compute centroid quadrature
			int = sum(self.elemAreas' .* F);

		end

		function int = threePointQuadrature(self,f)
		% Computes Gaussian three point quadrature on elements of the mesh.
		% Input f should be either a function_handle, a symbolic function, or a
		% vector of function evaluations on the midpoints of the edges of the
		% mesh. If:
		%
		%	(1) f is a function_handle, then f will be evaluated on the
		%		midpoints of the edges of the mesh.
		%	(2) f is a symbolic function, then f will be converted to a
		%		function_handle and evaluated on the midpoints of the edges of
		%		the mesh.
		%	(3) f is a vector, then f should be indexed so that f(i) stores the
		%		function evaluation on the i-th edge of the mesh, where the
		%		edges are ordered lexicographically. 
		%
		% Note: if f is a vector of function evaluations on the nodes of the
		% mesh (rather than on the midpoints), use nodalQuadrature.
		%
		% Note: Gaussian three point quadrature has order four convergence as
		% the mesh is refined. It is accurate to machine precision for
		% polynomials up to order 2. 

			% if f is a function
			if isa(f,'function_handle') || isa(f,'symfun')

				% convert symfun to function_handle
				if isa(f,'symfun'), f = matlabFunction(f); end

				% get edge midpoints
				edges = self.getMeshEdges;

				% get edge midpoints
				midpts = self.getMeshEdgeMidpoints(edges);

				% compute F on edge midpoints 
				F = f(midpts(:,1),midpts(:,2));

			% else if f is a vector of midpoint values
			elseif isvector(f)

				F = f; 

				% convert to column, if necessary
				if isrow(F), F = F'; end 

			end

			% compute quadrature
			elemAvg = sum([F(self.elemEdges(:,1)), ...
					F(self.elemEdges(:,2)),F(self.elemEdges(:,3))],2) / 3;
			int = sum(self.elemAreas' .* elemAvg);

		end


		%{
		% NOTE: I believe this function is mathematically equivalent to
		% nodalQuadrature, except it takes longer to run. I haven't done
		% testing to confirm this. But in the meantime, do not use this
		% function; rather, prefer nodalQuadrature.
		%----------------------------------------------------------------------%

		function int = threePointQuadrature_nodes(self,F_nodal)
		% Computes Gaussian three point quadrature on elements of the mesh.
		% Input F_nodal should be a vector of values computed on the nodes of
		% the mesh. It is assumed that F_nodal represents a solution to a
		% finite element problem using piecewise linear elements. Therefore,
		% the solution values on the midpoints of the triangles can be
		% calculated and used for Gaussian three point quadrature. 

			% get edge list
			edges = self.getMeshEdges;

			% compute F on midpoints
			F_midpts = (F_nodal(edges(:,1)) + F_nodal(edges(:,2))) / 2;

			% compute quadrature
			int = self.threePointQuadrature(F_midpts);
			
		end
		%----------------------------------------------------------------------%
		%}


		function IP = L2_IP_nodalQuadrature(self,arg1,arg2)

			% if arg1 is a function
			if isa(arg1,'function_handle') || isa(arg1,'symfun')

				% convert symfun to function_handle
				if isa(arg1,'symfun'), arg1 = matlabFunction(arg1); end

				% compute F on edge midpoints 
				Arg1 = arg1(self.Mesh.Nodes(1,:)',self.Mesh.Nodes(2,:)');

			% else if arg1 is a vector of nodal values
			elseif isvector(arg1)

				Arg1 = arg1; 

				% convert to column, if necessary
				if isrow(Arg1), Arg1 = Arg1'; end 

			end

			% if arg2 is a function
			if isa(arg2,'function_handle') || isa(arg2,'symfun')

				% convert symfun to function_handle
				if isa(arg2,'symfun'), arg2 = matlabFunction(arg2); end

				% compute F on edge midpoints 
				Arg2 = arg2(self.Mesh.Nodes(1,:)',self.Mesh.Nodes(2,:)');

			% else if arg1 is a vector of nodal values
			elseif isvector(arg2)

				Arg2 = arg2; 

				% convert to column, if necessary
				if isrow(Arg2), Arg2 = Arg2'; end 

			end

			% Compute the inner product
			IP = self.nodalQuadrature(Arg1 .* Arg2);

		end	

		function IP = L2_IP_threePointQuadrature(self,arg1,arg2)
		% Computes the L2 inner product using Gaussian three point quadrature.
		% Inputs arg1,arg2 should be either function_handles, symbolic functions, or
		% vectors of function evaluations on the midpoints of the edges of the
		% mesh. A mixture of these is also allowed. If:
		%
		%	(1) arg1 is a function_handle, then arg1 will be evaluated on the
		%		midpoints of the edges of the mesh.
		%	(2) arg1 is a symbolic function, then arg1 will be converted to a
		%		function_handle and evaluated on the midpoints of the edges of
		%		the mesh.
		%	(3) arg1 is a vector, then arg1 should be indexed so that arg1(i)
		%		stores the function evaluation on the i-th edge of the mesh,
		%		where the edges are ordered lexicographically. 
		%
		% The same is true for arg2.

			% if fxn-evals needed, get edge midpoints
			if isa(arg1,'function_handle') || ...
					isa(arg1,'symfun') || ...
					isa(arg2,'function_handle') || ...
					isa(arg2,'symfun')

				midpts = self.getMeshEdgeMidpoints;

			end

			% if arg1 is a function
			if isa(arg1,'function_handle') || isa(arg1,'symfun')

				% convert symfun to function_handle
				if isa(arg1,'symfun'), arg1 = matlabFunction(arg1); end

				% compute F on edge midpoints 
				Arg1 = arg1(midpts(:,1),midpts(:,2));

			% else if arg1 is a vector of midpoint values
			elseif isvector(arg1)

				Arg1 = arg1; 

				% convert to column, if necessary
				if isrow(Arg1), Arg1 = Arg1'; end 

			end

			% if arg2 is a function
			if isa(arg2,'function_handle') || isa(arg2,'symfun')

				% convert symfun to function_handle
				if isa(arg2,'symfun'), arg2 = matlabFunction(arg2); end

				% compute F on edge midpoints 
				Arg2 = arg2(midpts(:,1),midpts(:,2));

			% else if arg2 is a vector of midpoint values
			elseif isvector(arg2)

				Arg2 = arg2; 

				% convert to column, if necessary
				if isrow(Arg2), Arg2 = Arg2'; end 

			end

			% Compute the inner product
			IP = self.threePointQuadrature(Arg1 .* Arg2);

		end	

		function IP = L2_IP_threePointQuadrature_nodal(self,arg1,arg2)
		% Computes the L2 inner product using Gaussian three point quadrature.
		% Inputs arg1,arg2 should be either function_handles, symbolic functions, or
		% vectors of function evaluations on the nodes of the mesh.
		% A mixture of these is also allowed. 
		%
		% In particular, if arg1 is a vector of nodal values representing a
		% piecewise linear function, as would be the case for a Galerkin finite
		% element solution using linear Lagrange elements, then the values for
		% that solution on the edge midpoints can be exactly calculated, and
		% these can be passed to the usual threePointQuadrature function. The
		% purpose of threePointQuadrature_nodal is to convert a vectors of
		% nodal values to vectors of midpoint values, under this assumption. 
		%
		% WARNING: if arg1 is a vector of nodal values representing a function
		% that is *not* piecewise linear, then the midpoint values so computed
		% will not represent the true midpoint values, which may affect
		% accuracy or convergence of the quadrature. In this case, you should
		% use the L2_IP_ nodalQuadrature function instead, which expects
		% vectors of values computed on mesh nodes. 
		%
		% If:
		%
		%	(1) arg1 is a vector of nodal values, then these values will be
		%		averaged as appropriate to find the edge midpoint values. The
		%		vector of these midpoint values will be passed to
		%		L2_IP_threePointQuadrature function. 
		%	(2) arg1 is a function_handle or symbolic function, then arg1 will
		%		be passed directly to the L2_IP_threePointQuadrature function,
		%		which can handle function_handle or symbolic function inputs.
		%
		% All of the above statements are true for arg2 as well. 
		%
		% Note arg1 and arg2 can both be function_handles or symbolic
		% functions, in which case these arguments will be passed directly to
		% L2_IP_threePointQuadrature, essentially bypassing the current
		% function.

			% if finding the edge averages is needed, get edges
			if isa(arg1,'double') || isa(arg2,'double')
				edges = self.getMeshEdges;
			end

			if isa(arg1,'double')
				Arg1 = (arg1(edges(:,1)) + arg1(edges(:,2))) / 2;
			else
				Arg1 = arg1;
			end

			if isa(arg2,'double')
				Arg2 = (arg2(edges(:,1)) + arg2(edges(:,2))) / 2;
			else
				Arg2 = arg2;
			end

			IP = self.L2_IP_threePointQuadrature(Arg1,Arg2);

		end


		function norm = L2norm_threePointQuadrature(self,f)
		% Computes the L2 norm using Gaussian three point quadrature.
		% The norm is computed using the L2_IP_threePointQuadrature function.
		% Input f should be a function_handle, symbolic function, or vector of
		% function evaluations computed on the midpoints of the mesh.
		%
		% If:
		%
		%	(1) f is a function_handle or symbolic function, then f is passed
		%		directly to the L2_IP_threePointQuadrature function which can
		%		handle such inputs. 
		%	(2) f is a vector, then f should be indexed so that f(i) stores the
		%		function evaluation on the i-th edge of the mesh, where the
		%		edges are ordered lexicographically. 

			% compute the norm
			norm_squared = self.L2_IP_threePointQuadrature(f,f);
			norm = sqrt(norm_squared);

		end


		function norm = L2norm_threePointQuadrature_nodal(self,f)
		% Computes the L2 norm using Gaussian three point quadrature.
		% The norm is computed using the L2_IP_threePointQuadrature_nodal
		% function. Input f should be a function_handle, symbolic function, or
		% vector of function evaluations computed on the nodes of the mesh.
		%
		% In particular, if f is a vector of nodal values representing a
		% piecewise linear function, as would be the case for a Galerkin finite
		% element solution using linear Lagrange elements, then the values for
		% that solution on the edge midpoints can be exactly calculated. This
		% calculation is done by the L2_IP_threePointQuadrature_nodal function.
		% So the purpose of L2norm_threePointQuadrature_nodal is to route input
		% f to the corresponding L2_IP_threePointQuadrature_nodal function. 
		%
		% WARNING: if f is a vector of nodal values representing a function
		% that is *not* piecewise linear, then the midpoint values so computed
		% will not represent the true midpoint values, which may affect
		% accuracy or convergence of the quadrature. In this case, you should
		% use the L2norm_nodalQuadrature function instead, which expects
		% vectors of values computed on mesh nodes. 
		%
		% If:
		%
		%	(1) f is a function_handle or symbolic function, then f is passed
		%		directly to the L2_IP_threePointQuadrature_nodal function which
		%		can handle such inputs. 
		%	(2) f is a vector, then f should be indexed so that f(i) stores the
		%		function evaluation on the i-th node of the mesh, where the
		%		indexing is deteremined by self.Mesh.Nodes.

			% compute the norm
			norm_squared = self.L2_IP_threePointQuadrature_nodal(f,f);
			norm = sqrt(norm_squared);

		end

		% DEPRECATED! DO NOT USE IN THE FUTURE --------------------------------%
		function IP = L2_IP_piecewiseLinear(self,arg1,arg2)
		% arg1, arg2 should be a vectors of values on the nodes of the mesh,
		% representing a piecewise linear functions 

			% get mesh edges
			edges = self.getMeshEdges;

			% compute arg1,arg2 at midpoints (simple average)
			mdptVals1 = (arg1(edges(:,1)) + arg1(edges(:,2))) / 2;
			mdptVals2 = (arg2(edges(:,1)) + arg2(edges(:,2))) / 2;

			% compute L2 inner product using 3-pt gaussian quadrature
			IP = self.threePointQuadrature(mdptVals1 .* mdptVals2);

		end
		%----------------------------------------------------------------------%


		% DEPRECATED! DO NOT USE IN THE FUTURE --------------------------------%
		function int = L2norm_piecewiseLinear(self,F)
		% F should be a vector of values on the nodes of the mesh,
		% representing a piecewise linear functions 


			IP = L2_IP_piecewiseLinear(self,F,F);
			int = sqrt(IP);

		end
		%----------------------------------------------------------------------%


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

	end
end
