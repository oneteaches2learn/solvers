classdef NeumannBoundary < BoundaryCondition

	methods
		function obj = NeumannBoundary(condition,side)
			obj@BoundaryCondition(condition,side);
			[obj.valueLHS,obj.valueRHS] = obj.set_boundary_values(condition);
		end

		function [valueLHS, valueRHS] = set_boundary_values(obj, condition, transmissibility)
			if strcmp(obj.side,"L")
				valueRHS = condition;
			elseif strcmp(obj.side,"R")
				valueRHS = -condition;
			end

			valueLHS = 0;
		end

	end
end
