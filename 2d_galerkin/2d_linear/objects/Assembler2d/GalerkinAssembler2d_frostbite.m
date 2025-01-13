classdef GalerkinAssembler2d_frostbite < GalerkinAssembler2d_rxndiff

	properties
	end

	methods
		function self = GalerkinAssembler2d_frostbite
		end
	end

	methods (Static)
		function auxfun = assembleCoefficients(c,k,f,uInit)

			% set placeholder r
			r = 0;

			% call superclass method
			auxfun = assembleCoefficients@GalerkinAssembler2d_parabolic(c,k,r,f,uInit);
	
			%{
			% build r coefficient
			r = @(x1,x2,u,v) r_const * r_activ .* (u - v);

			% build dr_du coefficient (using product rule)
			r_activ_du = auxfun.cofs.r_activ_du;
			dr_du = @(x1,x2,u,v) r_const * (r_activ_du(u,v) .* (u - v) + r_activ(u,v));

			% store results
			auxfun.cofs.r = r;
			auxfun.cofs.dr_du = dr_du;
			%}

		end

		function auxfun = assembleReactionTerm(auxfun,r_const,r_activ,r_activ_du)
		% assembleReactionTerm assembles r and r_du from the comonents r_const,
		% r_activ, and possibly r_activ_du
		%
		% For the frostbite model, it is desireable to store r_const and r_activ
		% separately from r itself. In addition, r_active will frequently be an
		% anonymous function (not a symbolic function), so that r_activ_du
		% cannot be computed symbolically. In this case, r_activ_du must be
		% passed as an argument. If r_activ is a double or a symbolic function,
		% then r_activ_du will be computed symbolically.

			% define symbolic variables
			x = sym('x',[1 2],'real'); t = sym('t','real'); u = sym('u','real'); v = sym('v','real');

			% Process r_const (assuming r_const is a double)
			r_const = symfun(r_const,x);

			% Process r_activ
			if ~isa(r_activ,'function_handle')
				r_activ = symfun(r_activ,x);
			end

			% Process r_activ_du
			if isa(r_activ,'sym') || isa(r_activ,'symfun')
				r_activ_du = diff(r_activ,u);
			end

			% convert to matlab functions
			if isa(r_const,'symfun'), r_const = matlabFunction(r_const,"vars",[x t u v]); end		
			if isa(r_activ,'symfun'), r_activ = matlabFunction(r_activ,"vars",[x t u v]); end
			if isa(r_activ_du,'symfun'), r_activ_du = matlabFunction(r_activ_du,"vars",[x t u v]); end

			% Compute r
			r = @(x1,x2,t,u,v) r_const(x1,x2,t,u,v) .* r_activ(x1,x2,t,u,v) .* (u - v);

			% compute dr_du (using product rule)
			dr_du = @(x1,x2,t,u,v) r_const(x1,x2,t,u,v) .* (r_activ_du(x1,x2,t,u,v) .* (u - v) + r_activ(x1,x2,t,u,v));

			% store results
			auxfun.cofs.r = r;
			auxfun.cofs.dr_du = dr_du;
			auxfun.cofs.r_const = r_const;
			auxfun.cofs.r_activ = r_activ;
			auxfun.cofs.r_activ_du = r_activ_du;

		end
	end

end
