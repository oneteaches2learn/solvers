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

			% store variables
			elemNodes = self.Mesh.Elements';
			elemNodes(:,4) = elemNodes(:,1);

			% create lists of corresponding nodes
			ind = 1;
			for i = 1:self.nElem
				for j = 1:3
					
					% start nodes for each edge
					r(ind) = elemNodes(i,j);

					% end nodes for each edge
					s(ind) = elemNodes(i,j+1);

					ind = ind+1;

				end
			end

			% generate graph
			G = graph(r,s);
			
			% remove duplicate edges
			G = simplify(G);

		end			

		function edges = getMeshEdges(self)
		% get a list of all edges, ordered lexicographically

			% generate graph from mesh
			G = getMeshGraph(self);

			% store edges as an array
			edges = table2array(G.Edges);

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

			% store variables
			elemCoord = self.Mesh.Nodes';

			% compute midpoints
			startCoords = elemCoord(edges(:,1),:);
			endCoords = elemCoord(edges(:,2),:);
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

			pause();
			close();
		
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
		% Note: Gaussian three point quadrature has order two convergence as
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

		function int = L2norm_piecewiseLinear(self,F)
		% F should be a vector of values on the nodes of the mesh,
		% representing a piecewise linear functions 


			IP = L2_IP_piecewiseLinear(self,F,F);
			int = sqrt(IP);

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

	end
end
