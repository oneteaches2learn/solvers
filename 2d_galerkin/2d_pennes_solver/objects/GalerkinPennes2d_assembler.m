classdef GalerkinPennes2d_assembler < GalerkinParabolic2d_assembler

	properties
	end

	methods
		function self = GalerkinPennes2d_assembler
		end
	end

	methods (Static)
		function cofs = assembleCoefficients(p,k,r,uStar)

			% check function variables
			x = sym('x',[1 2],'real'); syms t;
			p = symfun(p,[x t]);
			k = symfun(k,[x t]);
			r = symfun(r,[x t]);
			uStar = symfun(uStar,[x t]);

			% convert to function_handles
			cofs.p = matlabFunction(p);
			cofs.k = matlabFunction(k);
			cofs.r = matlabFunction(r);
			cofs.uStar = matlabFunction(uStar);

		end
	end

	

end