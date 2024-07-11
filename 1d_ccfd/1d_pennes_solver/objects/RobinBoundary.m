classdef RobinBoundary < BoundaryCondition

	methods
		function obj = RobinBoundary(condition,side)
			obj@BoundaryCondition(condition,side);
			condition = obj.check_condition(condition);
			obj.condition = condition;
			[obj.valueLHS,obj.valueRHS] = obj.set_boundary_values(condition);
		end

		function condition = check_condition(obj,condition)
			err = ['For RobinBoundary, condition = [gamma,val] where q dot n = gamma(U - val)'];
			if ~isvector(condition) || length(condition) == 1
				error(err)
			end

			if iscolumn(condition)
				condition = condition';
			end

		end

		function [valueLHS, valueRHS] = set_boundary_values(obj, condition)
			valueLHS = condition(1);
			valueRHS = condition(1) * condition(2);
		end

	end
end
