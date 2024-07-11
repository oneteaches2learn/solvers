classdef Coefficient

	properties
		cellEdges
		cellCenters
		fun
	end

	methods
		function obj = Coefficient(domain,fun)
			obj.check_inputs(domain,fun);
			obj.fun = fun;
			obj.cellEdges = obj.evaluate_function(domain.cellEdges);
			obj.cellCenters = obj.evaluate_function(domain.cellCenters);
		end

		function check_inputs(obj,domain,fun)
			err = "Inputs for Coefficient must be a Domain1d object and a function handle";
			if ~isa(domain,'Domain') || ~isa(fun,'function_handle')
				error(err)
			end
		end

		function y = evaluate_function(obj,x)
			if nargin(obj.fun) == 0
				y = obj.fun() * ones(size(x));
			elseif nargin(obj.fun) == 1
				y = obj.fun(x);
			else
				err = "Check number of inputs for function_handle";
				error(err)
			end
		end

		function print(obj)
			cellEdges = full(obj.cellEdges')
			cellCenters = full(obj.cellCenters')
		end

	end
end



