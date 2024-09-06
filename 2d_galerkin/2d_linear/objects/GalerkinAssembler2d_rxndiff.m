classdef GalerkinAssembler2d_rxndiff < GalerkinAssembler2d_parabolic

	properties
	end

	methods
		function self = GalerkinAssembler2d_rxndiff
		end
	end

	methods (Static)
		function cofs = assembleCoefficients(p,k,r)

			% call superclass method
			cofs = assembleCoefficients@GalerkinAssembler2d_parabolic(p,k);

			% check function variables
			x = sym('x',[1 2],'real'); syms t;
			r = symfun(r,x);

			% convert to function_handles
			cofs.r = matlabFunction(r);

		end
	end

	

end
