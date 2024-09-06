classdef GalerkinAssembler2d_pennes < GalerkinAssembler2d_parabolic

	properties
	end

	methods
		function self = GalerkinAssembler2d_parabolic
		end
	end

	methods (Static)
		function cofs = assembleCoefficients(p,k,r,uStar)

			% call superclass method
			cofs = assembleCoefficients@GalerkinAssembler2d_parabolic(p,k);

			% check function variables
			x = sym('x',[1 2],'real'); syms t;
			r = symfun(r,x);
			uStar = symfun(uStar,x);

			% convert to function_handles
			cofs.r = matlabFunction(r);
			cofs.uStar = matlabFunction(uStar);

		end
	end

	

end
