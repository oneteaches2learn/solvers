classdef GalerkinPennes2d_assembler < GalerkinParabolic2d_assembler

	properties
	end

	methods
		function self = GalerkinPennes2d_assembler
		end
	end

	methods (Static)
		function cofs = assembleCoefficients(p,k,r,uStar)

			% call superclass method
			cofs = assembleCoefficients@GalerkinParabolic2d_assembler(p,k);

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
