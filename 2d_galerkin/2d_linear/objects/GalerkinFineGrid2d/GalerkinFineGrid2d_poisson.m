classdef GalerkinFineGrid2d_poisson < GalerkinFineGrid2d_elliptic

	properties
	end

	methods
		function self = GalerkinFineGrid2d_poisson(dom,auxfun,mmsparams,NameValueArgs)

			arguments
				dom 		Domain2d
				auxfun  	
				mmsparams 	MMSParams
				NameValueArgs.errType = "L2"
			end

			% call superclass constructor
			self@GalerkinFineGrid2d_elliptic(dom,auxfun,mmsparams,NameValueArgs)

		end

		function prob = solve(self,dom,cofs)

			prob = GalerkinSolver2d_poisson(dom,cofs);

		end

	end
end