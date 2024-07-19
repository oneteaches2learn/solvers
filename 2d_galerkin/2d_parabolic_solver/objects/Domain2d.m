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
