classdef GalerkinPoisson2d_mms < GalerkinParabolic2d_mms
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
		function self = GalerkinPoisson2d_mms(dom,auxfun,mmsparams,NameValueArgs)

			arguments
				dom 		Domain2d
				auxfun  	ManufacturedFunctions2d_poisson
				mmsparams 	MMSParams
				NameValueArgs.errType = "L2"
			end

			% create dummy TimeStepping object
			time = TimeStepping(0,0,0)


			% call superclass constructor
			self@GalerkinParabolic2d_mms(dom,time,auxfun,mmsparams,NameValueArgs)

		end

		function prob = solve(self,dom,time,cofs,uInit,f)

			prob = GalerkinPoisson2d_solver(dom,cofs,uInit,f);

		end

	end
end

