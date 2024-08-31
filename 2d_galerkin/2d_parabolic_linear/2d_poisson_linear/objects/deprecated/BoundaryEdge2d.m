classdef BoundaryEdge2d

	properties
		ID
		vertex1
		vertex2
		outwardNormal
		boundaryType
		boundaryCondition
		nodes
		nNodes
	end

	methods
		function self = BoundaryEdge2d(vertex1,vertex2,n)

			self.vertex1 = vertex1;
			self.vertex2 = vertex2;
			self.outwardNormal = n;

		end

		function self = setNodes()
		end

		function self = setBCs()
		end

	end
end
