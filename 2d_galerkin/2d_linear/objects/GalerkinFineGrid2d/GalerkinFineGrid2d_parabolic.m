classdef GalerkinFineGrid2d_parabolic < GalerkinFineGrid2d
% GALERKINFINEGRID2D_PARABOLIC(DOM,AUXFUN,MMSPARAMS,NAMEVALUEARGS) runs an MMS
% test on a parabolic PDE
%
% author: Tyler Fara				Date: July 11, 2024
%-----------------------------------------------------------------------------%
%
% required arguments
%	dom			Domain2d object
%	auxfun		
%	mmsparams	MMSParams object
%
% optional arguments
%	errType
%
% properties
%	NEEDS TO BE FILLED IN
%
%-----------------------------------------------------------------------------%

	properties
	end

	methods
		function self = GalerkinFineGrid2d_parabolic(dom,auxfun,mmsparams,NameValueArgs)

			% call superclass constructor
			self@GalerkinFineGrid2d(dom,auxfun,mmsparams,NameValueArgs);

		end

	end
end