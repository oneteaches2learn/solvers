classdef NewtonGalerkinMMS2d_ginzburglandau < GalerkinMMS2d_elliptic
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
		function self = NewtonGalerkinMMS2d_ginzburglandau(dom,auxfun,mmsparams,NameValueArgs)

			arguments
				dom 		Domain2d
				auxfun  	ManufacturedFunctions2d_ginzburglandau
				mmsparams 	MMSParams
				NameValueArgs.errType = "L2"
			end

			% call superclass constructor
			self@GalerkinMMS2d_elliptic(dom,auxfun,mmsparams,NameValueArgs)

		end

		function prob = solve(self,dom,cofs)

			%prob = GalerkinSolverNR2d_poisson(dom,cofs);
			prob = NewtonGalerkinSolver2d_ginzburglandau(dom,cofs);

		end

	end
end
