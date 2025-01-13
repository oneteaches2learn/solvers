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

