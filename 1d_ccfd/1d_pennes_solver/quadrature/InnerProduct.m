classdef InnerProduct

	properties
		quadrature
	end

	methods
		function obj = InnerProduct(quad)
			obj.check_inputs(quad);
			obj.quadrature = quad;
		end

		function val = compute(obj,f1,f2,x)
			err = "Inputs must be vectors of same length";
			if length(f1) ~= length(f2), error(err), end

			product = f1 .* f2;
			val = obj.quadrature.compute(x,product);
		end

		function check_inputs(obj,quad)
			err = "Input must be a QuadratureRule object.";
			if ~isa(quad,'QuadratureRule')
				error(err)
			end
		end
	end
end
