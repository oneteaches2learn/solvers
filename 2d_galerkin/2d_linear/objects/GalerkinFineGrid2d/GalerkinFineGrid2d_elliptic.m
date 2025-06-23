classdef GalerkinFineGrid2d_elliptic < GalerkinFineGrid2d
% GALERKINFINEGRID2D_ELLIPTIC(DOM,AUXFUN,MMSPARAMS,NAMEVALUEARGS) runs an MMS
% test on a parabolic PDE
%
% author: Tyler Fara				Date: May 20, 2025
%-----------------------------------------------------------------------------%
%
% required arguments
%	dom			Domain2d object
%	auxfun		AuxFunctions2d object
%	mmsparams	MMSParams object
%
% optional arguments
%	NameValueArgs
%
% properties
%	NEEDS TO BE FILLED IN
%
%-----------------------------------------------------------------------------%

	methods
		function self = GalerkinFineGrid2d_elliptic(dom,auxfun,mmsparams,NameValueArgs)

			% call superclass constructor
			self@GalerkinFineGrid2d(dom,auxfun,mmsparams,NameValueArgs);

		end

	end
end