classdef GalerkinRxndiff2d_assembler < GalerkinParabolic2d_assembler

	properties
	end

	methods
		function self = GalerkinRxndiff2d_assembler
		end
	end

	methods (Static)
		function cofs = assembleCoefficients(p,k,r)

			% call superclass method
			cofs = assembleCoefficients@GalerkinParabolic2d_assembler(p,k);

			% check function variables
			x = sym('x',[1 2],'real'); syms t;
			r = symfun(r,x);

			% convert to function_handles
			cofs.r = matlabFunction(r);

		end
	end

	

end
