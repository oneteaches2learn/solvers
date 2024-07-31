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

			% generate graph from mesh
			G = getMeshGraph(self);

			% store edges as an array
			edges = table2array(G.Edges);

		end

		function midptCoords = getMeshEdgeMidpoints(self,edges)

			% store variables
			elemCoord = self.Mesh.Nodes';

			% compute midpoints
			startCoords = elemCoord(edges(:,1),:);
			endCoords = elemCoord(edges(:,2),:);
			midptCoords = (startCoords + endCoords) / 2;

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
			
			% if f is function_handle
			if isa(f,'function_handle')

				F = f(self.Mesh.Nodes(1,:)',self.Mesh.Nodes(2,:)');

			% if f is symbolic function
			elseif isa(f,'sym')

				f = matlabFunction(f);
				F = f(self.Mesh.Nodes(1,:)',self.Mesh.Nodes(2,:)');

			% if f is vector of nodal values
			elseif isvector(f)

				if isrow(f), F = f'; else, F = f; end

			elseif isdouble(f)

				F = f * ones(self.nNodes,1);

			end

			% loop over elements
			int = 0;
			elemNodes = self.Mesh.Elements';
			for i = 1:self.nElem

				% nodal quadrature on element
				int = int + self.elemAreas(i) * sum(F(elemNodes(i,:))) / 3;

			end
		end

		function int = centroidQuadrature(self,f)
			
			% compute element centroids
			centroids = zeros(self.nElem,2);
			elemNodes = self.Mesh.Elements';
			elemCoord = self.Mesh.Nodes';
			for i = 1:self.nElem
				centroids(i,:) = sum(elemCoord(elemNodes(i,:),:)) / 3;
			end

			% if f is function_handle
			if isa(f,'function_handle')

				F = f(centroids(:,1),centroids(:,2));

			% if f is symbolic function
			elseif isa(f,'sym')

				f = matlabFunction(f);
				F = f(centroids(:,1),centroids(:,2));

			% if f is vector of nodal values
			elseif isvector(f)

				if isrow(f), F = f'; else, F = f; end

			elseif isdouble(f)

				F = f * ones(self.nElem,1);

			end

			% loop over elements
			int = 0;
			elemInd = self.Mesh.Elements';
			for i = 1:self.nElem

				% nodal quadrature on element
				int = int + self.elemAreas(i) * F(i);

			end
		end

		function int = threePointQuadrature(self,f)

			edges = self.getMeshEdges;

			if isa(f,'function_handle') || isa(f,'symfun')

				% convert symfun to function_handle
				if isa(f,'symfun'), f = matlabFunction(f); end

				% get edge midpoints
				midpts = self.getMeshEdgeMidpoints(edges);

				% compute F on edge midpoints 
				F = f(midpts(:,1),midpts(:,2));

			elseif isvector(f)

				F = f; 

				% convert to column, if necessary
				if isrow(F), F = F'; end 

			end

			% store F as sparse, symmetric, weighted adjacency matrix
			F = sparse([edges(:,1);edges(:,2)], ...
					[edges(:,2);edges(:,1)], ...
					[F(:);F(:)]);

			% loop on elements
			elemNodes = self.Mesh.Elements';
			elemNodes(:,4) = elemNodes(:,1);
			int = 0;
			for i = 1:self.nElem

				% loop on edges
				tot = 0;
				for j = 1:3
					tot = tot + F(elemNodes(i,j),elemNodes(i,j+1));
				end
				
				% compute three pt quadrature on element
				int = int + self.elemAreas(i) * tot/3;

			end

		end

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
