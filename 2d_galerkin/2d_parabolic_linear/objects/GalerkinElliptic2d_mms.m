classdef GalerkinElliptic2d_mms < Galerkin2d_mms
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
		function self = GalerkinElliptic2d_mms(dom,auxfun,mmsparams,NameValueArgs)

			% call superclass constructor
			self@Galerkin2d_mms(dom,auxfun,mmsparams,NameValueArgs);

		end

	end
end

