classdef Domain1d

	properties
		cellEdges
		cellWidths
		cellCenters
		cellNum
		length
	end

	methods
		function obj = Domain1d(cellEdges)
			cellEdges = check_input(obj,cellEdges);

			obj.cellEdges = cellEdges;
			obj.cellWidths = diff(obj.cellEdges);
			obj.cellNum = length(obj.cellEdges) - 1;
			obj.cellCenters = get_cell_centers(obj);
			obj.length = sum(obj.cellWidths);
		end

		function cellEdges = check_input(obj,cellEdges)
			[n,m] = size(cellEdges);

			err = "Input for Domain1d needs to be a vector of cell edges";
			if ~isvector(cellEdges), error(err), end

			if n ~= 1
				cellEdges = cellEdges';
			end
		end

		function cellCenters = get_cell_centers(obj)
			N = obj.cellNum;
			dx = obj.cellWidths;
			cellCenters = zeros(1,N);
			cellCenters = obj.cellEdges(1:N) + 0.5 * dx;
		end

	end
end
