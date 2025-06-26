classdef GalerkinAssembler2d_pennes < GalerkinAssembler2d_parabolic

	properties
	end

	methods
		function self = GalerkinAssembler2d_parabolic
		end
	end

	methods (Static)
		function auxfun = assembleCoefficients(c,k,r,uStar,f,uInit)

			% call superclass method
			auxfun = assembleCoefficients@GalerkinAssembler2d_parabolic(c,k,r,f,uInit);

			% store blood temperature
			x = sym('x',[1 2],'real');

			if ~isa(uStar,'function_handle')
				auxfun.cofs.uStar = matlabFunction(symfun(uStar,x));
			else
				auxfun.cofs.uStar = uStar;
			end

		end
	end

	

end
