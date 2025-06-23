classdef GalerkinFineGrid2d_heat < GalerkinFineGrid2d_parabolic

	properties
	end

	methods
		function self = GalerkinFineGrid2d_heat(dom,auxfun,mmsparams,NameValueArgs)

			arguments
				dom 		Domain2d
				auxfun  	
				mmsparams 	MMSParams
				NameValueArgs.errType = "Linfty(L2)"
			end

			% call superclass constructor
			self@GalerkinFineGrid2d_parabolic(dom,auxfun,mmsparams,NameValueArgs)

		end

		function prob = solve(self,dom,cofs)

			prob = GalerkinSolver2d_heat(dom,cofs);

		end

	end
end
