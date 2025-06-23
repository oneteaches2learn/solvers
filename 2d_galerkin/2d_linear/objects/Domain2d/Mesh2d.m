classdef Mesh2d
% Mesh2d is an interface with MATLAB's FEMesh object

	properties
		base
		p
		h
		nNodes
		nElems
		nEdges
		nodes
		elements
		edges
		edgeLengths
		effectiveNodes
		effectiveElems
		unusedNodes
		unusedElems
		elementCoords
		elementEdges
		midpoints
		centroids
		areas
		quality
	end

	properties (Hidden)
		Mesh
	end

	methods
		% CONSTRUCTOR
		function self = Mesh2d(dl,p,base)

			% construct mesh from dl, p, and base
			if nargin == 3

				% store inputs
				self.p = p;
				self.base = base;
				self.h = base^-p;

				% generate the mesh
				self.Mesh = self.setMesh(dl,p,base);

				% precompute areas to save computation time later
				self.areas = self.get_areas;
			end

			% create mesh from FEMesh object
			if nargin == 1 && isa(dl,'pde.FEMesh')
				FEmsh = dl;
				self.Mesh = FEmsh;
				self.base = 2;
				self.p = 1;
				self.h = 1;

				% precompute areas to save computation time later
				self.areas = self.get_areas;
			end

		end


		% SETTERS
		function mesh = setMesh(self,dl,p,base)

			% create mesh
			geo = fegeometry(dl);
			geo = geo.generateMesh(Hmax=base^-p,GeometricOrder='linear');
			mesh = geo.Mesh;

		end

		function self = setEffectiveNodes(self)
			
			% store node data
			self.effectiveNodes = [1:1:self.nNodes];
			self.unusedNodes = [];

		end

		function self = setEffectiveElements(self)

			% store element data
			self.effectiveElems = [1:1:self.nElems];
			self.unusedElems = [];

		end

		function self = setFreeNodes(self)

			% store free nodes
			self.freeNodes = setdiff(1:self.mesh.nNodes,self.boundaryNodes.D);

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

		function vals = get.edgeLengths(self)

			% get edge data
			edg = self.edges;
			nod = self.nodes;

			% get edge coordinates
			x1 = nod(edg(:,1),1);
			y1 = nod(edg(:,1),2);
			x2 = nod(edg(:,2),1);
			y2 = nod(edg(:,2),2);

			% compute edge lengths
			vals = sqrt((x1 - x2).^2 + (y1 - y2).^2);
			vals = vals(:);

		end

		function val = get.h(self)

			% get mesh size
			val = max(self.edgeLengths);

		end

		function vals = get.elementCoords(self)
		% returns an nElem x 3 x 2 array of element coordinates
		%    the first index is the element number, the second index is the node
		%    number, the third index is the x or y coordinate.
		

			% get node data
			elemNodes = self.elements;
			elemCoord = self.nodes;

			% isolate coordanates of each node for each element
			nodes1 = elemCoord(elemNodes(:,1),:);
			nodes2 = elemCoord(elemNodes(:,2),:);
			nodes3 = elemCoord(elemNodes(:,3),:);

			% concatenate node coordinates for each element
			vals = cat(1,nodes1,nodes2,nodes3);
			vals = reshape(vals,[],3,2);

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

		function vals = get.quality(self)
		% The metric of quality used here is min(2 * r / R), where r is the
		% in-radius and R is the circum-radius of the triangle. The value of 1
		% indicates an equilateral triangle, while values less than 1 indicate a
		% degenerate triangle.

			% get element data
			elemNodes = self.elements;
			elemCoord = self.nodes;

			% isolate coordanates of each node for each element
			nodes1 = elemCoord(elemNodes(:,1),:);
			nodes2 = elemCoord(elemNodes(:,2),:);
			nodes3 = elemCoord(elemNodes(:,3),:);

			% get edge lengths
			edge1 = sqrt((nodes1(:,1) - nodes2(:,1)).^2 + ...
							(nodes1(:,2) - nodes2(:,2)).^2);
			edge2 = sqrt((nodes1(:,1) - nodes3(:,1)).^2 + ...
							(nodes1(:,2) - nodes3(:,2)).^2);
			edge3 = sqrt((nodes2(:,1) - nodes3(:,1)).^2 + ...
							(nodes2(:,2) - nodes3(:,2)).^2);

			% get in-radius
			r = 2 * self.areas ./ (edge1 + edge2 + edge3);

			% get circum-radius
			R = (edge1 .* edge2 .* edge3) ./ ...
				(4 * self.areas);

			% compute quality of each element
			vals = min(2 * r ./ R);

		end

		function vals = get_areas(self)

			% if mesh was generated by MATLAB's PDE toolbox
			if isa(self.Mesh,'pde.FEMesh')
				[temp,vals] = area(self.Mesh);
				vals = vals';

			% else mesh was generated by some other means
			else
        		% Extract node indices for each element
            	elem_nodes = self.elements;  % nElems x 3 matrix
            
				% Retrieve the coordinates of each node for all elements
				x1 = self.nodes(elem_nodes(:,1), 1);
				y1 = self.nodes(elem_nodes(:,1), 2);
				
				x2 = self.nodes(elem_nodes(:,2), 1);
				y2 = self.nodes(elem_nodes(:,2), 2);
				
				x3 = self.nodes(elem_nodes(:,3), 1);
				y3 = self.nodes(elem_nodes(:,3), 2);
				
				% Compute the area using the shoelace formula
				% Area = |(x1(y2 - y3) + x2(y3 - y1) + x3(y1 - y2))| / 2
				areas = abs(x1 .* (y2 - y3) + ...
							x2 .* (y3 - y1) + ...
							x3 .* (y1 - y2)) / 2;
				
				% Ensure that areas is a column vector
				vals = areas(:);
			end

			% note: if self.Mesh is an FEmesh object, then computation of the
			% areas is about 10x faster. 
			
		end

		function angles = computeTriangleAngles(self)
		% computeTriangleAngles  Compute interior angles (in degrees) of all triangles.
		%
		%   angles = computeTriangleAngles(nodes, elements)
		%
		%   Inputs:
		%     nodes    – an (N_nodes × 2) array of node coordinates:
		%                nodes(i,:) = [ x_i,  y_i ]
		%     elements – an (N_elems × 3) array of node indices (1-based) for each triangle:
		%                elements(j,:) = [ i1,  i2,  i3 ]
		%
		%   Output:
		%     angles   – an (N_elems × 3) array of interior angles (in degrees).  
		%                For each triangle j, 
		%                   angles(j,1) = angle at node i1,
		%                   angles(j,2) = angle at node i2,
		%                   angles(j,3) = angle at node i3.
		%
		%   Example:
		%     nodes = [ 0, 0;
		%               1, 0;
		%               0.5, 0.5;
		%               0, 1 ];
		%     elements = [ 1, 2, 3;
		%                  1, 3, 4 ];
		%     ang = computeTriangleAngles(nodes, elements)
		%     % returns:
		%     %   ang =
		%     %     45    45    90
		%     %     45    90    45

			% Extract nodes and elements from the mesh
			nodes = self.nodes;
			elements = self.elements;
			
			% Number of elements
			nElems = size(elements, 1);

			% Pre‐allocate arrays
			angle1 = zeros(nElems,1);
			angle2 = zeros(nElems,1);
			angle3 = zeros(nElems,1);

			% Extract coordinates of each vertex in each triangle
			p1 = nodes(elements(:,1), :);   % (nElems × 2)
			p2 = nodes(elements(:,2), :);   % (nElems × 2)
			p3 = nodes(elements(:,3), :);   % (nElems × 2)

			% Compute edge vectors for each triangle
			v12 = p2 - p1;   % vector from vertex 1 to vertex 2
			v13 = p3 - p1;   % vector from vertex 1 to vertex 3
			v23 = p3 - p2;   % vector from vertex 2 to vertex 3

			v21 = -v12;      % vector from vertex 2 to vertex 1
			v31 = -v13;      % vector from vertex 3 to vertex 1
			v32 = -v23;      % vector from vertex 3 to vertex 2

			% Compute lengths of each edge
			l12 = sqrt(sum(v12.^2, 2));   % length of edge 1–2
			l13 = sqrt(sum(v13.^2, 2));   % length of edge 1–3
			l23 = sqrt(sum(v23.^2, 2));   % length of edge 2–3

			% Avoid division by zero by enforcing a small positive minimum
			eps_val = 1e-15;
			l12 = max(l12, eps_val);
			l13 = max(l13, eps_val);
			l23 = max(l23, eps_val);

			% --- Angle at vertex 1: between edges v12 and v13 ---
			% dot(v12, v13) / (|v12| * |v13|)
			cosAngle1 = ( v12(:,1).*v13(:,1) + v12(:,2).*v13(:,2) ) ./ (l12 .* l13);
			cosAngle1 = min(max(cosAngle1, -1.0), 1.0);  % clip to [−1, 1]
			angle1 = acosd(cosAngle1);

			% --- Angle at vertex 2: between edges v21 and v23 ---
			cosAngle2 = ( v21(:,1).*v23(:,1) + v21(:,2).*v23(:,2) ) ./ (l12 .* l23);
			cosAngle2 = min(max(cosAngle2, -1.0), 1.0);
			angle2 = acosd(cosAngle2);

			% --- Angle at vertex 3: between edges v31 and v32 ---
			cosAngle3 = ( v31(:,1).*v32(:,1) + v31(:,2).*v32(:,2) ) ./ (l13 .* l23);
			cosAngle3 = min(max(cosAngle3, -1.0), 1.0);
			angle3 = acosd(cosAngle3);

			% Stack the three angles into an (nElems × 3) array
			angles = [ angle1,  angle2,  angle3 ];
		end

		% SETTERS
		function self = set.nodes(self,nodes)

			self.Mesh.Nodes = nodes';

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

		function h = plotSolution(self,sol)
		% plots the solution on the mesh
		%
		% inputs
		% sol - solution vector, nNodes x 1
		%
		% outputs
		% h - figure handle
		%
			% check if solution is a column vector
			if size(sol,2) > 1
				error("Solution must be a column vector")
			end

			% plot solution
			h = surf(self.nodes(1,:),self.nodes(2,:),sol);
			%h = pdeplot(self.Mesh,'XYData',sol,'ZData',sol);
			%title('Solution on Mesh');
			%xlabel('X-axis');
			%ylabel('Y-axis');
			colorbar;
			
		end

	end

end


