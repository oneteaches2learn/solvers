classdef TrapezoidalRule < QuadratureRule
%TRAPEZOIDALRULE represents the trapezoidal rule for integration of function f.
%	TRAPEZOIDALRULE inherits from QuadratureRule. The OBJ.COMPUTE(X,F) method
%	computes the trapezoidal rule for F, a vector of function evaluations, on
%	domain X, which is a vector of grid points. TRAPEZOIDALRULE expects that F =
%	f(X).

	methods
		function val = compute(obj,x,f)
			obj.check_inputs(x,f);

			dx = diff(x);
			M  = length(dx);
			val = 0; 

			for i = 1:M
				val = val + 0.5 * (f(i) + f(i+1)) * dx(i);
			end
		end

		function check_inputs(obj,f,x)
			err = "Inputs must be two vectors of same length.";
			if ~isvector(x) || ~isvector(f)
				error(err)
			end
			if length(x) ~= length(f)
				error(err)
			end
		end
	end
end

