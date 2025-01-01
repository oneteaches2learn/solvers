classdef CoupledNewtonMMS2d_parabolic < CoupledNewtonMMS2d
% GALERKINELLIPTIC2D_MMS(DOM,TIME,AUXFUN,MMSPARAMS,NAMEVALUEARGS) runs an MMS
% test on a parabolic PDE
%
% author: Tyler Fara				Date: July 11, 2024
%-----------------------------------------------------------------------------%
%
% required arguments
%	dom			Domain2d object
%	time		TimeStepping object
%	auxfun		AuxFunctions2d object
%	mmsparams	MMSParams object
%
% optional arguments
%	errType
%
% properties
%	NEEDS TO BE FILLED IN
%
%-----------------------------------------------------------------------------%

	methods
		function self = CoupledNewtonMMS2d_parabolic(dom,auxfun,mmsparams,NameValueArgs)

			% call superclass constructor
			self@CoupledNewtonMMS2d(dom,auxfun,mmsparams,NameValueArgs);

		end

		function problems = solveManufacturedProblems(self,varargin)

			% unpack variables
			cofs = self.auxFunctions.functionHandles;
			base = self.mmsParams.base;
			pmin = self.mmsParams.pmin;
			pmax = self.mmsParams.pmax;
			tOff = self.mmsParams.timeOffset;
			tFac = self.mmsParams.timeFactor;
			region = self.mmsParams.effectiveRegion;

			% run mms test
			fprintf('MMS Test Begun\n')
			fprintf('Solving Problems\n')
			ind = 1;
			for p = pmin:pmax
				
				% successively refine mese
				fprintf(' p = %i solved:',p); tic;

				dom_p = self.configureDomain(p); 
				prob_p = self.solve(dom_p,cofs);

				% store results
				if self.mmsParams.demo == 0, problems{ind} = prob_p;
				else problems{1} = prob_p;
				end

				% prepare next trial
				ind = ind + 1;
				executionTime = toc;
				fprintf(' %f s\n',executionTime)

			end
		end


	end
end

