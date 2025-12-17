classdef Coefficients

	methods
		function self = Coefficients
		end
	end


	% COEFFICIENT BANK
	methods (Static)
		function val = ramp_activation(u,v,cofs)

			arguments
				u
				v
				cofs.v_min= 0;
				cofs.v_max = 1;
				cofs.u_min = 0;
				cofs.u_max = 1;
				cofs.gamma = 1;
			end
			
			% define individual ramps
			ramp_v = @(v) Coefficients.ramp(v,lowerBound=cofs.v_min,upperBound=cofs.v_max);
			ramp_u = @(u) Coefficients.ramp(u,lowerBound=cofs.u_min,upperBound=cofs.u_max);

			% combine maps
			cof = @(u,v) (ramp_u(u) .* ramp_v(v).^cofs.gamma);
			%cof = @(u,v) (ramp_u(u) .* ramp_v(v).^cofs.gamma + 0.1);
			val = cof(u,v);

		end	

		function val = ramp_activation_du(u,v,cofs)

			arguments
				u
				v
				cofs.v_min = 0;
				cofs.v_max = 1;
				cofs.u_min = 0;
				cofs.u_max = 1;
				cofs.gamma = 1;
			end

			% define individual functions
			ramp_v = @(v) Coefficients.ramp(v,lowerBound=cofs.v_min,upperBound=cofs.v_max);
			ramp_u_du = @(u) Coefficients.ramp_du(u,lowerBound=cofs.u_min,upperBound=cofs.u_max);

			% combine maps
			cof = @(u,v) ramp_u_du(u) .* ramp_v(v).^cofs.gamma;
			val = cof(u,v);

		end	

		function val = logistic_activation(u,v,cofs)

			arguments
				u
				v
				cofs.u_L = 1;
				cofs.u_k = 1;
				cofs.u_0 = 0;
				cofs.v_L = 1;
				cofs.v_k = 1;
				cofs.v_0 = 0;
				cofs.gamma = 1;
			end

			% define individual functions
			logistic_v = @(v) Coefficients.logistic(v,L=cofs.v_L,k=cofs.v_k,u_0=cofs.v_0);
			logistic_u = @(u) Coefficients.logistic(u,L=cofs.u_L,k=cofs.u_k,u_0=cofs.u_0);

			% combine maps
			cof = @(u,v) logistic_u(u) .* logistic_v(v).^cofs.gamma;
			val = cof(u,v);

		end

		function val = logistic_activation_du(u,v,cofs)
			
			arguments
				u
				v
				cofs.u_L = 1; % supremum of logistic function
				cofs.u_k = 1; % steepness of logistic function
				cofs.u_0 = 0; % midpoint of logistic function
				cofs.v_L = 1;
				cofs.v_k = 1;
				cofs.v_0 = 0;
				cofs.gamma = 1;
			end

			% define individual functions
			logistic_v = @(v) Coefficients.logistic(v,L=cofs.v_L,k=cofs.v_k,u_0=cofs.v_0);
			logistic_u_du = @(u) Coefficients.logistic_du(u,L=cofs.u_L,k=cofs.u_k,u_0=cofs.u_0);

			% combine maps
			cof = @(u,v) logistic_u_du(u) .* logistic_v(v).^cofs.gamma;
			val = cof(u,v);

		end

		function [r,s] = coupling_coefficients(type,data)
		% coupling_coefficients(type,data) returns the coupling coefficients
		% r and s for the given type and data.
		%
		% INPUTS
		%	type - string, type of coupling coefficients, e.g. 'R1', 'R2', 'R3'
		%	data - struct, data for the coupling coefficients
		%
		% OUTPUTS
		%	r - struct, coupling coefficients for the reaction term
		%	s - struct, coupling coefficients for the source term
		%
		% The coupling coefficients are used to compute the coefficients for the
		% coupling terms between the PDE and ODE. These coupling terms are of
		% the form:
		% 
		% 		r(u,v) = r_const * r_activ(u,v) * (u - v)   
		% 		s(u,v) = s_const * s_activ(u,v) * (v - u)
		%
		% Coefficients for the PDE are assembled by GalerkinAssembler2d, which
		% expects double r_const and function_handles r_activ and r_activ_du
		% (the derivative of r_active with respect to u). The ODE object expects
		% a single function_handle combining s_const and s_activ. So,
		% COUPLING_COEFFICIENTS assembles these objects and stores them as
		% fields of r and s. 
		

			% STORE COEFFICIENTS
			if strcmp(type,'R1') || strcmp(type,'R2') || strcmp(type,'R3')

				% store ramp data
    			v_min = data.v_min; 
				v_max = data.v_max; 
				u_min = data.u_min; 
				u_max = data.u_max; 
				gamma = data.gamma;
				r_const = data.r_const;
				s_const = data.s_const;

			elseif strcmp(type,'S1') || strcmp(type,'S2') || strcmp(type,'S3')

				% store logistic data
				u_L = data.u_L; 
				u_k = data.u_k; 
				u_0 = data.u_0;
				v_L = data.v_L; 
				v_k = data.v_k; 
				v_0 = data.v_0;
				gamma = data.gamma; 
				r_const = data.r_const;
				s_const = data.s_const;

			elseif strcmp(type,'Stolwijk')
			% note: right now, you have no way of passing other parameters to
			% the Stolwijk model. You should rethink how you're doing this. 
			
				% store Stolwijk data
				r_const = data.r_const;
				s_const = data.s_const;

			else
				error('Unknown type of coupling coefficients: %s',type);

			end

			% RAMP ACTIVATION COEFFICIENTS
			% if cof depends on v only
			if strcmp(type,'R1')

				% compute r_activ
				r_activ = @(x1,x2,t,u,v) Coefficients.ramp(v, ...
											lowerBound=v_min, upperBound=v_max);
				
				% compute r_activ_du
				r_activ_du = @(x1,x2,t,u,v) 0;

				% compute s_activ
				s_activ = @(u,v) r_activ(0,0,0,u,v);

			% if cof depends on u only
			elseif strcmp(type,'R2')

				% compute r_activ
				r_activ = @(x1,x2,t,u,v) Coefficients.ramp(u, ...
											lowerBound=u_min, upperBound=u_max);
				
				% compute r_activ_du
				r_activ_du = @(x1,x2,t,u,v) Coefficients.ramp_du(u, ...
                                lowerBound=u_min, upperBound=u_max);

				% compute s_activ
				s_activ = @(u,v) r_activ(0,0,0,u,v);

			% if cof depends on u and v
			elseif strcmp(type,'R3')

				% compute r_activ
				r_activ = @(x1,x2,t,u,v) Coefficients.ramp_activation(u,v, ...
					v_min=v_min,v_max=v_max,u_min=u_min,u_max=u_max,gamma=gamma);

				% compute r_activ_du
				r_activ_du = @(x1,x2,t,u,v) Coefficients.ramp_activation_du(u,v, ...
					v_min=v_min,v_max=v_max,u_min=u_min,u_max=u_max,gamma=gamma);
				
				% compute s_activ
				s_activ = @(u,v) r_activ(0,0,0,u,v);

			end


			% LOGISTIC ACTIVATION COEFFICIENTS
			% if cof depends on v only
			if strcmp(type,'S1')

				% compute r_activ
				r_activ = @(x1,x2,t,u,v) Coefficients.logistic(v, ...
											L=v_L,k=v_k,u_0=v_0);
				
				% compute r_activ_du
				r_activ_du = @(x1,x2,t,u,v) 0;

				% compute s_activ
				s_activ = @(u,v) r_activ(0,0,0,u,v);

			% if cof depends on u only
			elseif strcmp(type,'S2')

				% compute r_activ
				r_activ = @(x1,x2,t,u,v) Coefficients.logistic(u, ...
											L=u_L,k=u_k,u_0=u_0);
				
				% compute r_activ_du
				r_activ_du = @(x1,x2,t,u,v) Coefficients.logistic_du(u, ...
								L=u_L,k=u_k,u_0=u_0);

				% compute s_activ
				s_activ = @(u,v) r_activ(0,0,0,u,v);

			% if cof depends on u and v
			elseif strcmp(type,'S3')

				% compute r_activ
				r_activ = @(x1,x2,t,u,v) Coefficients.logistic_activation(u,v, ...
					u_L=u_L,u_k=u_k,u_0=u_0,v_L=v_L,v_k=v_k,v_0=v_0,gamma=gamma);

				% compute r_activ_du
				r_activ_du = @(x1,x2,t,u,v) Coefficients.logistic_activation_du(u,v, ...
					u_L=u_L,u_k=u_k,u_0=u_0,v_L=v_L,v_k=v_k,v_0=v_0,gamma=gamma);
				
				% compute s_activ
				s_activ = @(u,v) r_activ(0,0,0,u,v);

			end


			% STOLWIJK MODEL
			if strcmp(type,'Stolwijk')

				% compute r_activ
				r_activ = @(x1,x2,t,u,v) Coefficients.omega_s(u,v);

				% compute r_activ_du
				r_activ_du = @(x1,x2,t,u,v) Coefficients.omega_s_du(u,v);

				% compute s_activ
				s_activ = @(u,v) r_activ(0,0,0,u,v);

			end

			% PACKAGE RESULTS
			r.r_const = r_const;
			r.r_activ = r_activ;
			r.r_activ_du = r_activ_du;
			s.s_const = data.s_const;
			s.s_activ = s_activ;
			s.s_func  = @(u,v) s_const * r_activ(0,0,0,u,v);

		end
	end


	% BUILDING BLOCKS
	methods (Static)
		function val = rLU(u)
		% rLU(u) returns val = u for u >= 0, or val = 0 otherwise.
		%
		% NOTE: rLU = rectified linear unit

			val = max(u,0);

		end

		function val = ramp(u,cofs)
		% ramp(u) returns val = u / (upperBound - lowerBound) for lowerBound <= u
		% 	<= upperBound, or val = 0 for u < lowerBound, or val = 1 for u >
		% 	upperBound.
		% 
		% lowerBound and upperBound are passed as name-value pairs and have
		% default values lowerBound = 0 and upperBound = 1.  
		% 
		% syntax:
		% 	y = ramp_function(u)
		% 	y = ramp_function(u,lowerBound=-1)
		% 	y = ramp_function(u,upperBound=6)
		% 	y = ramp_function(u,lowerBound=-1,upperBound=6)

			arguments
				u
				cofs.lowerBound = 0;
				cofs.upperBound = 1;
			end

			rampLow  = @(u) Coefficients.rLU(u - cofs.lowerBound);
			rampHigh = @(u) Coefficients.rLU(u - cofs.upperBound);
			ramp = @(u) (rampLow(u) - rampHigh(u)) / (cofs.upperBound - cofs.lowerBound);

			val = ramp(u);

		end

		function val = ramp_du(u,cofs)
		% d_ramp(u) returns val = 1 / (upperBound - lowerBound) for lowerBound <= u
		% 	<= upperBound, and val = 0 otherwise. It is the derivative of ramp(u).
		%
		% lowerBound and upperBound are passed as name-value pairs and have
		% default values lowerBound = 0 and upperBound = 1.

			arguments
				u
				cofs.lowerBound = 0;
				cofs.upperBound = 1;
			end
			
			val = double(u >= cofs.lowerBound & u <= cofs.upperBound) / (cofs.upperBound - cofs.lowerBound);

		end

		function val = box(u,cofs)
		% box(x) returns val = 1 for lowerBound <= u <= upperBound, and val = 0 else.
		%
		% lowerBound and upperBound are passed as name-value pairs and have
		% default values lowerBound = 0 and upperBound = Inf.

			arguments
				u
				cofs.lowerBound = 0;
				cofs.upperBound = Inf;
			end

			val = double(u >= coefs.lowerBound & u <= coefs.upperBound);

		end

		function val = heaviside(u)
		% heaviside(u) returns val = 1 for u >= 0 and val = 0 otherwise.

			val = double(u >= 0);

		end
		
		function val = logistic(u,cofs)
		% logistic(u,cofs) returns val = L/(1 + exp(-k * (u - u_0)).
		%
		% Coefficients L, k, and u_0 are passed as name-value pairs and have
		% default values L = 1, k = 1, u_0 = 0.
		%
		% syntax:
		% 	y = logistic_function(u)
		% 	y = logistic_function(u,u_0=2)
		% 	y = logistic_function(u,L=2,k=2,u_0=2)

			arguments
				u
				cofs.L = 1;
				cofs.k = 1;
				cofs.u_0 = 0;
			end
		
			val = cofs.L ./ (1 + exp(-cofs.k * (u - cofs.u_0)));

		end
		
		function val = logistic_du(u,cofs)
		% logistic_du(u,cofs) returns val = L*k*exp(-k*(u - u_0))/(1 + exp(-k*(u - u_0))^2.
		%
		% Coefficients L, k, and u_0 are passed as name-value pairs and have
		% default values L = 1, k = 1, u_0 = 0.
		%
		% syntax:
		% 	y = logistic_function(u)
		% 	y = logistic_function(u,u_0=2)
		% 	y = logistic_function(u,L=2,k=2,u_0=2)

			arguments
				u
				cofs.L = 1;
				cofs.k = 1;
				cofs.u_0 = 0;
			end
		
			% explicit formulation of derivative
			%val = cofs.L * cofs.k * exp(-cofs.k * (u - cofs.u_0)) ./ (1 + exp(-cofs.k * (u - cofs.u_0)).^2);

			% formulation using logistic function
			f = @(u) Coefficients.logistic(u,L=cofs.L,k=cofs.k,u_0=cofs.u_0);
			val = f(u) .* (1 - f(u) / cofs.L) * cofs.k;

		end
	end

	% STOLWIJK MODEL
	methods (Static)

		function val = e_hy(v,e_hySET) 

			arguments
				v
				e_hySET = 37;
			end

			val = v - e_hySET;
		end

		function val = e_s(u,e_sSET)
			
			arguments
				u
				e_sSET = 33;
			end

			val = u - e_sSET;		

		end

		function val = DI(e_hy,e_s,C1,C2)
			
			arguments
				e_hy
				e_s
				C1 = 28424;
				C2 = 4870;
			end

			%val = Coefficients.rLU(C1 * e_hy + C2 * e_s);
			val = 1e-6 * Coefficients.rLU(C1 * e_hy + C2 * e_s);

		end

		function val = DI_du(e_hy,e_s,C1,C2)
			
			arguments
				e_hy
				e_s
				C1 = 28424;
				C2 = 4870;
			end

			%val = C2 * Coefficients.heaviside(C1 * e_hy + C2 * e_s);
			val = 1e-6 * C2 * Coefficients.heaviside(C1 * e_hy + C2 * e_s);

		end


		function val = CS(e_hy,e_s,C3,C4)
			
			arguments
				e_hy
				e_s
				C3 = 1.1;
				C4 = 3.3;
			end

			%val = Coefficients.rLU(-C3 * e_hy - C4 * e_s);
			val = 1e-6 * Coefficients.rLU(-C3 * e_hy - C4 * e_s);

		end

		function val = CS_du(e_hy,e_s,C3,C4)
			
			arguments
				e_hy
				e_s
				C3 = 1.1;
				C4 = 3.3;
			end

			%val = -C4 * Coefficients.heaviside(-C3 * e_hy - C4 * e_s);
			val = -1e-6 * C4 * Coefficients.heaviside(-C3 * e_hy - C4 * e_s);

		end

		function val = omega_s(u,v,omega_0,alpha_DI,cofs)

			arguments
				u
				v
				omega_0  = 4.39e-4; 
				alpha_DI = 0.134;
				cofs.C1  = 28424;
				cofs.C2  = 4870;
				cofs.C3  = 1.1;
				cofs.C4  = 3.3;
				cofs.e_hySET = 37;
				cofs.e_sSET = 34;
			end

			e_hy = Coefficients.e_hy(v,cofs.e_hySET);
			e_s  = Coefficients.e_s(u,cofs.e_sSET);
			DI   = Coefficients.DI(e_hy,e_s,cofs.C1,cofs.C2);
			CS   = Coefficients.CS(e_hy,e_s,cofs.C3,cofs.C4);

			val = (omega_0 + alpha_DI .* DI) ./ (1 + CS);

		end

		function val = omega_s_du(u,v,omega_0,alpha_DI,cofs)

			arguments
				u
				v
				omega_0  = 4.39e-4;
				alpha_DI = 0.134;
				cofs.C1  = 28424;
				cofs.C2  = 4870;
				cofs.C3  = 1.1;
				cofs.C4  = 3.3;
				cofs.e_hySET = 37;
				cofs.e_sSET = 34;
			end

			e_hy  = Coefficients.e_hy(v,cofs.e_hySET);
			e_s   = Coefficients.e_s(u,cofs.e_sSET);
			DI    = Coefficients.DI(e_hy,e_s,cofs.C1,cofs.C2);
			DI_du = Coefficients.DI_du(e_hy,e_s,cofs.C1,cofs.C2);
			CS    = Coefficients.CS(e_hy,e_s,cofs.C3,cofs.C4);
			CS_du = Coefficients.CS_du(e_hy,e_s,cofs.C3,cofs.C4);

			val = (alpha_DI .* DI_du .* (1 + CS) - (omega_0 + alpha_DI .* DI) .* CS_du) ./ (1 + CS).^2;

		end
	end



	% VARIABLE CHECKING / SETTING METHODS
	methods (Static)

		function [f,t,U,V] = checkVariables(f)

			% store function type data
			isCoupled     = num2str(Coefficients.isCoupled(f));
			isNonlinear   = num2str(Coefficients.isNonlinear(f));
			isTimeVarying = num2str(Coefficients.isTimeVarying(f));
			code = strcat(isCoupled,isNonlinear,isTimeVarying);
			code = str2num(code);

			% Update function signature and store variables and/or set dummy variables
			switch code
				% Spatially varying only
				case 000
					f = @(x1,x2,t,u,v)(f(x1,x2));

				% Time varying
				case 001
					f = @(x1,x2,t,u,v)(f(x1,x2,t));
				
				% Nonlinear
				case 010
					f = @(x1,x2,t,u,v)(f(x1,x2,u));
				
				% Time varying and nonlinear
				case 011
					f = @(x1,x2,t,u,v)(f(x1,x2,t,u));
				
				% Coupled
				case 100
					f = @(x1,x2,t,u,v)(f(x1,x2,v));

				% Coupled and time varying
				case 101
					f = @(x1,x2,t,u,v)(f(x1,x2,t,v));

				% Coupled and nonlinear
				case 110
					f = @(x1,x2,t,u,v)(f(x1,x2,u,v));

				% Coupled, time varying, and nonlinear
				case 111
					f = @(x1,x2,t,u,v)(f(x1,x2,t,u,v));
			end

		end

		%{
		function [f,t,u] = checkVariables(f)
		% checkVariables(f) converts f to an anonymous function handle in the
		% variables x, t, and u.
		% 
		% INPUTS
		%	f - double, sym or symfun, or function handle
		%
		% OUTPUTS
		%	f - function handle, f(x1,x2,t,u)
		%	t - double, time variable
		%	u - double, nonlinear variable
		%
		% In this package, the most general need is for f to be a nonlinear
		% function of space and time, where space is two dimensional. To support
		% this, all function handles are presumed to be functions of x, t, and
		% u. This assumption is made even when solving linear problems, or
		% stationary problems. In the case that f is constant, linear, or
		% stationary, the additional variables are added and set to zero. In the
		% case that f is not an anonymous function handle, i.e. if f is a
		% double, sym, or symfun, then f is first converted to a function
		% handle.

			if ~Coefficients.isNonlinear(f) && ~Coefficients.isTimeVarying(f)
				f = @(x1,x2,t,u)(f(x1,x2));
				t = 0;
				u = 0;
			elseif ~Coefficients.isNonlinear(f) && Coefficients.isTimeVarying(f)
				f = @(x1,x2,t,u)(f(x1,x2,t));
				u = 0;
			elseif Coefficients.isNonlinear(f) && ~Coefficients.isTimeVarying(f)
				f = @(x1,x2,t,u)(f(x1,x2,u));
				t = 0;
			end
		end
		%}

		function f = wrap2d(f)
		% wrap2d(f) ensures that f returns an output of the same size as the
		% input x1.
		%
		% INPUTS
		%	f - function handle, f(x1,x2,t,u)
		%
		% OUTPUTS
		%	f - function handle, f(x1,x2,t,u) + zeros(size(x1))
		%
		% In this package, f will usually be used to return a function at each
		% node of some mesh; i.e. the output of f will usually be a matrix.
		% However, if f is constant in the variable x, then MATLAB's
		% broadcasting rules will cause f to return a scalar. To ensure that f
		% returns a matrix, wrap2d adds an appropriately sized matrix of zeros
		% to the output of f.
		
			f = @(x1,x2,t,u)(f(x1,x2,t,u) + zeros(size(x1)));
		end


		function result = isTimeVarying(f)
		% isTimeVarying(f) returns true if the function handle f is time varying.
		%
		% INPUTS
		%	f - double, sym or symfun, or function handle
		%
		% OUTPUTS
		%	result - logical, true if f is time varying, false otherwise
		%
		% isTimeVarying checks if input f is a function of variable t. To do so,
		% isTimeVarying must be converted into an anonymous function handle. If f
		% is a double, then f represents a constant function and is not time
		% varying. In this case, f is converted to a symfun in the spatial
		% variable x. If f is a sym or symfun, then f is converted to a function
		% handle. The function handle is then checked for the presence of the
		% variable t.

			if isa(f,'double')
				x = sym('x',[1 2]);
				f = symfun(f,x);
			end
			
			if isa(f,'sym') || isa(f,'symfun')
				x = sym('x',[1 2]);
				f = matlabFunction(f,x);
			end

			result = Coefficients.hasVariable(f,'t');
		end

		function result = isNonlinear(f)
		% isNonlinear(f) returns true if the function handle f is nonlinear in
		% the variable u.
		%
		% INPUTS
		%	f - double, sym or symfun, or function handle
		%
		% OUTPUTS
		%	result - logical, true if f is nonlinear in u, false otherwise
		%
		% isNonlinear checks if input f is a function of variable u. To do so,
		% isNonlinear must be converted into an anonymous function handle. If f
		% is a double, then f represents a constant function and is nonlinear in
		% u. In this case, f is converted to a symfun in the spatial variable x.
		% If f is a sym or symfun, then f is converted to a function handle. The
		% function handle is then checked for the presence of the variable u.
		
			if isa(f,'double')
				x = sym('x',[1 2]);
				f = symfun(f,x);
			end

			if isa(f,'sym') || isa(f,'symfun')
				f = matlabFunction(f);
			end
			
			result = Coefficients.hasVariable(f,'u');
		end

		function result = isCoupled(f)
		% isNonlinear(f) returns true if the function handle f is a function of
		% the ODE variable v.
		%
		% INPUTS
		%	f - double, sym or symfun, or function handle
		%
		% OUTPUTS
		%	result - logical, true if f is a function of v, false otherwise
		%
		% isNonlinear checks if input f is a function of variable v. To do so,
		% isNonlinear must be converted into an anonymous function handle. If f
		% is a double, then f represents a constant function and is not a
		% function of v. In this case, f is converted to a symfun in the spatial
		% variable x.  If f is a sym or symfun, then f is converted to a
		% function handle. The function handle is then checked for the presence
		% of the variable v.
		
			if isa(f,'double')
				x = sym('x',[1 2]);
				f = symfun(f,x);
			end

			if isa(f,'sym') || isa(f,'symfun')
				f = matlabFunction(f);
			end
			
			result = Coefficients.hasVariable(f,'v');
		end

		function result = hasVariable(f,checkVar)
		% hasVariable(f,checkVar) returns true if the function handle f is a
		% function of the variable checkVar.
		%
		% INPUTS
		%	f - double, sym or symfun, or function handle
		%	checkVar - string, variable to check for
		%
		% OUTPUTS
		%	result - logical, true if f is a function of checkVar, false otherwise
		%
		% In this package, functions are presumed to be functions of x, t, and
		% u. checkVar could, presumably, be some other variable. But hasVariable
		% has been designed and tested for checkVar = 't' or 'u'.
		
			if isa(f,'double')
				x = sym('x',[1 2]);
				f = symfun(f,x);
			end
			
			if isa(f,'sym') || isa(f,'symfun')
				f = matlabFunction(f,x);
			end
			
			usedVars = Coefficients.getFunctionVariables(f);
			result = any(strcmp(usedVars,checkVar));
		end

		function vars = getFunctionVariables(fh)
		% getFunctionVariables(fh) returns the variables of function handle fh
		%
		% INPUTS
		%	fh - function handle
		%
		% OUTPUTS
		%	vars - cell array of strings, variables of fh
		%
		% getFunctionVariables uses regex to collect variable names from
		% function_handle fh. For example, if fh is a function of x, t, and u,
		% then vars will be {'x1', 'x2', 't','u'}. Presumably,
		% getFunctionVariables could check for variables other than x, t, and u.
		% But getFunctionVariables has been designed and tested for x, t, and u.
		
			% use regex to collect variable names from function_handle
			code = func2str(fh);
			pattern = '@\((.*?)\)';
			tokens = regexp(code,pattern,'tokens');

			if ~isempty(tokens)
				vars = strsplit(tokens{1}{1},',');
				vars = strtrim(vars);

			else
				vars = {};

			end
		end
	end

end

