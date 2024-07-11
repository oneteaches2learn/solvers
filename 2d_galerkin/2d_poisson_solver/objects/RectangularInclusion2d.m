classdef RectangularInclusion2d

	properties
		vertices
		xBounds
		yBounds
		edges
	end

	methods
		function self = RectangularInclusion2d(gd_col)

			% set vertices
			self.vertices = self.geometryDescriptionColumnToVertices(gd_col);

			% set x- and y-bounds
			self.xBounds = [self.vertices(1,1) self.vertices(2,1)];
			self.yBounds = [self.vertices(1,2) self.vertices(3,2)];

			% set edges
			self.edges = self.setEdges();

		end

		function edges = setEdges(self)

			% set up vertex matrix for cycle
			vert = self.vertices;
			vert = [vert; vert(1,:)];

			% set up outward normal
			n = [0 1; -1 0; 0 -1; 1 0];

			for i = 1:4
				edges{i} = BoundaryEdge2d(vert(i,:),vert(i+1,:),n(i,:));
			end

		end

		function self = setBoundaryTypes(self,boundaryTypes)

			for i = 1:4
				edges{i}.boundaryType = boundaryTypes{i};
			end

		end

		function self = setBoundaryConditions(self,boundaryConditions)

			for i = 1:4
				edges{i}.boundaryCondition = boundaryConditions{i};
			end

		end

		function vert = geometryDescriptionColumnToVertices(self,gd_col)

			vert_x = gd_col(3:6);
			vert_y = gd_col(7:10);

			vert = [vert_x vert_y];

		end
	end

end
