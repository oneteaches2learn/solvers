classdef NewtonGalerkinMMS2d < GalerkinMMS2d
% GALERKIN2D_MMS(DOM,AUXFUN,MMSPARAMS,NAMEVALUEARGS) runs an MMS
% test on a parabolic PDE
%
% author: Tyler Fara				Date: October 12, 2024
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

	properties
	end

	methods
		% CONSTRUCTOR 
		function self = NewtonGalerkinMMS2d(dom,auxfun,mmsparams,NameValueArgs)

			% call superclass constructor
			self@GalerkinMMS2d(dom,auxfun,mmsparams,NameValueArgs);
		
		end


		% UTILITY FUNCTIONS
		function dom = manufactureBoundaryConditions(self,dom,auxfun,varargin);

			% unpack coefficients
			uTrue = self.auxFunctions.uTrue;
			if isa(self,'GalerkinMMS2d_parabolic')
				u_t = self.auxFunctions.u_t;
			end
			q = self.auxFunctions.q;
			nEdges = dom.boundary.nEdges;

			% set symbolic variables
			if isa(self,'GalerkinMMS2d_parabolic')
				x = sym('x',[1 2]); syms t;
				vars = [x t];
			else
				x = sym('x',[1 2]);
				vars = x;
			end

			% TEMPORARY: VARIABLE HANDLING, ALWAYS INCLUDES U...SHOULD IT STAY THIS WAY??
			x = sym('x',[1 2]); syms t; syms u;
			vars = [x u];

			% set edge normal vectors
			dom = self.setEdgeNormalVectors_outerBoundary(dom);
			if isa(dom,'Domain2d_punctured')
				dom = self.setEdgeNormalVectors_inclusions(dom);
			end

			% Manufacture Dirichlet BC
			u_d = uTrue;
			u_d = matlabFunction(u_d);

			% assign boundary functions
			for i = 1:nEdges

				% assign Dirichlet BC
				if dom.boundary.edges(i).boundaryType == 'D'
					dom.boundary.edges(i).boundaryCondition = u_d;

				% assign Neumann BC
				elseif dom.boundary.edges(i).boundaryType == 'N'

					% store coefficients
					n_i = dom.boundary.edges(i).outwardNormal;
					u_N = symfun(auxfun.u_N,vars);

					% compute residual
					n_i = symfun(n_i,vars);
					q_i = symfun(sum(q.*n_i),vars);

					% attempt 1: doesn't work, compose is too hard to work with
					%s_i = q_i - compose(u_N,uTrue,u,x);

					% attempt 2: works, but isn't extensible to time-varying problem
					s_i = q_i - u_N(x(1),x(2),uTrue(x(1),x(2)));

					% assemble BC
					g_i  = u_N + s_i;
					dg_i = diff(g_i,u);
					g_i  = matlabFunction(g_i);
					dg_i = matlabFunction(dg_i);

					% store result
					dom.boundary.edges(i).boundaryCondition = g_i;
					dom.boundary.edges(i).boundaryCondition_ddu = dg_i;

				% assign Robin BC
				elseif dom.boundary.edges(i).boundaryType == 'R'
					alpha_i = symfun(1.0,vars);
					n_i = symfun(dom.boundary.edges(i).outwardNormal,vars);
					g_i = symfun(uTrue - sum(q .* n_i) / alpha_i,vars);
					alpha_i = matlabFunction(alpha_i);
					g_i = matlabFunction(g_i);
					dom.boundary.edges(i).boundaryCondition = {alpha_i,g_i};

				% assign dynamic BC
				elseif dom.boundary.edges(i).boundaryType == 'T'
					alpha_i = symfun(1.0,vars);
					beta_i  = symfun(1.0,vars);
					gamma_i = symfun(1.0,vars);
					n_i = symfun(dom.boundary.edges(i).outwardNormal,vars);
					g_i = symfun((alpha_i * uTrue + beta_i * u_t - sum(q .* n_i)) / gamma_i,vars);
					alpha_i = matlabFunction(alpha_i);
					beta_i  = matlabFunction(beta_i);
					gamma_i = matlabFunction(gamma_i);
					g_i = matlabFunction(g_i);
					dom.boundary.edges(i).boundaryCondition = {alpha_i,beta_i,gamma_i,g_i};
				end
			end
		end

	end
end

