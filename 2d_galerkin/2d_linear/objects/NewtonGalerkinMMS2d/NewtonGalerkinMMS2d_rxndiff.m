classdef NewtonGalerkinMMS2d_rxndiff < NewtonGalerkinMMS2d_parabolic
% DOCUMENTATION NEEDED! In the meantime, see the documentation for
% poissonFD1d_finegrid, which is a very similar object. 
%
% Notes
%	(1) Use the member funtion self.instantiateFineGridTest to create a fine
%	grid test using the same input data as the corresponding mms test. This can
%	be used to compare results for an mms versus fine grid test.
%-----------------------------------------------------------------------------%

	properties
	end

	methods
		function self = NewtonGalerkinMMS2d_rxndiff(dom,auxfun,mmsparams,NameValueArgs)

			arguments
				dom 		Domain2d
				auxfun  	ManufacturedFunctions2d_rxndiff
				mmsparams 	MMSParams
				NameValueArgs.errType = "Linfty(L2)"
			end

			% call superclass constructor
			self@NewtonGalerkinMMS2d_parabolic(dom,auxfun,mmsparams,NameValueArgs)

		end

		function prob = solve(self,dom,cofs)

			prob = NewtonGalerkinSolver2d_rxndiff(dom,cofs);

		end

	end
end
