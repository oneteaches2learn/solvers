classdef Domain

	properties
		cellEdges
		cellWidths
		cellCenters
		cellNum
		length
	end

	methods
		function obj = Domain(cellEdges)
			cellEdges = check_input(obj,cellEdges);

			obj.cellEdges = cellEdges;
			obj.cellWidths = diff(obj.cellEdges);
			obj.cellNum = length(obj.cellEdges) - 1;
			obj.cellCenters = get_cell_centers(obj);
			obj.length = sum(obj.cellWidths);
		end

		function cellEdges = check_input(obj,cellEdges)
			[n,m] = size(cellEdges);

			err = "Input for Domain needs to be a vector of cell edges";
			if ~isvector(cellEdges), error(err), end

			if isrow(cellEdges), cellEdges = cellEdges'; end
		end

		function cellCenters = get_cell_centers(obj)
			N = obj.cellNum;
			dx = obj.cellWidths;
			cellCenters = zeros(1,N);
			cellCenters = obj.cellEdges(1:N) + 0.5 * dx;
		end

		function print(obj)
			cellEdges = obj.cellEdges'
			cellCenters = obj.cellCenters'
			cellWidths = obj.cellWidths'
			cellNum = obj.cellNum
			length = obj.length
		end
	end
end
