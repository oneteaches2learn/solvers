classdef Mesh2d
% Mesh2d is an interface with MATLAB's FEMesh object

	properties
		nNodes
		nElems
		nEdges
		nodes
		elements
		edges
		elementEdges
		midpoints
		centroids
		areas
	end

	properties (Hidden)
		Mesh
	end

	methods
		% CONSTRUCTOR
		function self = Mesh2d(mesh)

			self.Mesh = mesh;
			%{
			self.nNodes = self.nNodes;
			self.nElems = self.nElems;

			self.nEdges = self.nEdges;
			self.nodes = self.nodes;
			self.elements = self.elements;
			self.edges = self.edges;
			self.elementEdges = self.elementEdges;
			self.midpoints = self.midpoints;
			self.centroids = self.centroids;
			%}
			self.areas = self.get_areas;
		end

		% GETTERS
		function nodes = get.nodes(self)

			nodes = self.Mesh.Nodes';

		end

		function elems = get.elements(self)

			elems = self.Mesh.Elements';

		end

		function edges = get.edges(self)

			% get list of all element edges
			elemNodesOrig = self.elements';
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

		function num = get.nNodes(self)

			num = size(self.nodes,1);
			
		end

		function num = get.nElems(self)

			num = size(self.elements,1);

		end

		function num = get.nEdges(self)

			num = size(self.edges,1);
			
		end

		function vals = get.elementEdges(self)

			% get edge and element data
			edges     = self.edges;
			elemNodes = sort(self.elements,2);
			
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

			vals = elemEdges;

		end

		function midpts = get.midpoints(self)
			
			midpts = (self.nodes(self.edges(:,1),:) + ...
								self.nodes(self.edges(:,2),:)) / 2;

		end

		function vals = get.centroids(self)

			% get centroid coordinates
			elemNodes = self.elements;
			elemCoord = self.nodes;

			% isolate coordanates of each node for each element
			nodes1 = elemCoord(elemNodes(:,1),:);
			nodes2 = elemCoord(elemNodes(:,2),:);
			nodes3 = elemCoord(elemNodes(:,3),:);

			% compute average of element nodes for each element
			vals = (nodes1 + nodes2 + nodes3) / 3;

		end

		function vals = get_areas(self)

			[temp,vals] = area(self.Mesh);
			vals = vals';

		end


		% PLOTTERS
		function h = plot(self,NameValueArgs)

			arguments
				self
				NameValueArgs.NodeLabels string = "off"
				NameValueArgs.NodeFontSize double = 18
				NameValueArgs.ElementLabels string = "off"
			end
			x = NameValueArgs;

			% Plot mesh
			hold on
			h = pdemesh(self.Mesh,ElementLabels=x.ElementLabels);

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
			
			h = gcf;

		end

	end

end


