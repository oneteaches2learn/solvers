classdef MMSParams
% MMSPARAMS(BASE,NAMEVALUEARGS) contains parameters for parabolic MMS test
%
% author: Tyler Fara						date: July 11, 2024
%-----------------------------------------------------------------------------%
%
% required inputs
% 	base		integer, base for grid size, i.e. h = base^p
%
% optional inputs
%	demo		integer, default 0, means run full MMS test
%	pmin		integer, minimum exponent for MMS test
%	pmax		integer, maximum exponent for MMS test
%	timeOffset	integer, offset for time grid 
%	timeFactor	integer, scaling factor for time grid
%
%-----------------------------------------------------------------------------%
%
% notes
%	(1) The MMS test will construct SPATIAL domains with mesh-refinement 
%			
%				h = base^pmin    through    h = base^max.
%
%	By default, pmin = 1 and pmax = 4. These values can be overwritten by
%	passing pmin and pmax as name-value pairs. 
%
%	(2) The MMS test will construct TIME domains with mesh-refinement
%
%				dt = base^(timeFactor * (p - timeOffset)).
%
%	To obtain optimal convergence rates, it is usually necessary that the
%	time-grid be refined twice as fast as the space-grid. Hence, timeFactor =
%	2, by default. This default can be overwritten by passing the desired
%	timeFactor as a name-value pair. The default timeFactor (i.e. timeFactor =
%	2) may result in many time steps, causing error computation to take a very
%	long time, so the timeOffset allows to reduce the number of time steps,
%	while still scaling the number of time steps in a way that results in
%	optimal convergence. By default, timeOffset = 1. The timeOffset can be
%	overwritten by passing timeOffset as a name-value pair. 
%
%	(3) Demo-mode: Often it is desirable to run only one trial of the MMS test.
%	(For example, when testing the solver, one might desire to use the MMS test
%	to manufacture boundary conditions, RHS, and coarse mesh on some
%	domain...but not to run successively refined grids, compute error, etc.
%	Hence the "demo" property. By default, demo = 0, which runs an MMS test
%	from h = base^pmin through h = base^pmax. However, setting demo to any
%	other positive integer (done by passing demo as a name-value pair) will put
%	the MMS test into demo-mode, wherein pmin = pmax = demo; i.e. the value
%	passed for demo will be used for the power, and only one trial will be run
%	with h = base^demo as the mesh refinement.
%
%-----------------------------------------------------------------------------%
%
% FUNCTIONS
% [pmin,pmax] = setOrderBounds(self,NameValueArgs)
% 	Sets pmin and pmax. By default, pmin = 1 and pmax = 4. 
%
%	If user instantiates MMSParams with other pmin or pmax values (passed as
%	name-value pairs to the constructor) then setOrderBounds will use these
%	values instead. Alternatively, user may instantiate MMSParams with demo=p
%	as a name-value pair (where p is a positive integer). In this case,
%	setOrderBounds will set pmin = pmax = p, causing the MMS test to run just
%	one trial with spatial refinement h = base^p and time refinement dt =
%	timeFactor * base^(p - timeOffset).
%
%-----------------------------------------------------------------------------%

	properties
		base
		demo
		pmin
		pmax
		timeOffset
		timeFactor
		effectiveRegion
		meshInclusions
		quadType
	end

	properties (Dependent)
		nTrials
	end

	methods
		function self = MMSParams(base,NameValueArgs)

			arguments
				base
				NameValueArgs.demo = 0
				NameValueArgs.pmin = 1
				NameValueArgs.pmax = 4
				NameValueArgs.timeOffset = 1
				NameValueArgs.timeFactor = 2
				NameValueArgs.effectiveRegion = 'Omega_eps'
				NameValueArgs.meshInclusions = 'off'
				NameValueArgs.quadType = 'threePoint'
			end

			% set properties
			self.base = base;
			self.demo = NameValueArgs.demo;
			self.timeOffset = NameValueArgs.timeOffset;
			self.timeFactor = NameValueArgs.timeFactor;
			self.effectiveRegion = NameValueArgs.effectiveRegion;
			self.meshInclusions = NameValueArgs.meshInclusions;
			self.quadType = NameValueArgs.quadType;
			[self.pmin,self.pmax] = self.setOrderBounds(NameValueArgs);
			
		end

		function [pmin,pmax] = setOrderBounds(self,NameValueArgs)

			% if not in demo-mode, set use default or user-passed values
			if self.demo == 0
				pmin = NameValueArgs.pmin;
				pmax = NameValueArgs.pmax;

			% else, if in demo-mode, use user-passed demo value
			else
				pmin = self.demo;
				pmax = self.demo;
			end

		end

		function nTrials = get.nTrials(self)
			nTrials = self.pmax - self.pmin + 1;
		end

	end

end
