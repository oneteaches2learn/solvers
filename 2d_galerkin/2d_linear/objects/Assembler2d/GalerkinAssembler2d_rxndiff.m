classdef GalerkinAssembler2d_rxndiff < GalerkinAssembler2d_parabolic

	properties
	end

	methods
		function self = GalerkinAssembler2d_rxndiff
		end
	end

	methods (Static)
		function auxfun = assembleCoefficients(c,k,r,f,uInit)

			% call superclass method
			auxfun = assembleCoefficients@GalerkinAssembler2d_parabolic(c,k,r,f,uInit);

		end
	end

	

end
