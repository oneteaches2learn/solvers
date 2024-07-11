classdef DirichletBoundary < BoundaryCondition

	methods
		function obj = DirichletBoundary(condition,side,transmissibility)
			obj@BoundaryCondition(condition,side);

			if ~isa(transmissibility,'Transmissibility') 
				error("DirichletBoundary(condition,side,transmissibility) requires Transmissibility object.")
			end

			[obj.valueLHS,obj.valueRHS] = obj.set_boundary_values(condition, transmissibility);
		end

		function [valueLHS, valueRHS] = set_boundary_values(obj, condition, transmissibility)
			if strcmp(obj.side,"L")
				trans = transmissibility.cellEdges(1);
			else
				n = length(transmissibility.cellEdges);
				trans = transmissibility.cellEdges(n);
			end

			valueLHS = trans;
			valueRHS = trans * condition;
		end

	end
end
