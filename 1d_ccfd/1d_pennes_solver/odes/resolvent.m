function v = resolvent(v_in,dt,c_B,vLower,vUpper)
%RESOLVENT(V_IN,DT,C_B,VLOWER,VUPPER) computes the resolvent for the ODE
%
%			v' + c_B(v - S(u)) + lambda = mu
%
%	where lambda is a Lagrange multiplier enforcing some constraint.
%
% Author: Tyler Fara			Date: April 3, 2024
%-----------------------------------------------------------------------------%
% Inputs
%	v_in	double or vector of doubles
%	dt		double, represents time step size
%	c_B		double, nonnegative
%	vLower	double or NaN, lower constraint, i.e. v >= vLower
%	vUpper  double or NaN, upper constraint, i.e. v <= vUpper
%
% Outputs
%	v		double or vector of doubles
%-----------------------------------------------------------------------------%
% Notes
%	(1) RESOLVENT is a general function that can be turned into a specific,
%	anonymous function R(v_in) by choosing DT, C_B, VLOWER, VUPPER in advance.
%	Supposing these variables have been assigned in advance, then the syntax is
%
%		R = @(v)(resolvent(v,dt,c_B,vLower,vUpper);
%
%	In this way, you can create custom resolvents enforcing lower, upper, or
%	both kinds of constraints. And, these custom anonymous functions can be
%	passed as arguments to other functions.
%-----------------------------------------------------------------------------%

	denom = 1 + dt * c_B;

	% CASE 1: with upper constraint, no lower constraint
	if isnan(vLower) && ~isnan(vUpper)
		for i = 1:length(v_in)
			if v_in(i) < vUpper* denom
				v(i) = v_in(i) / denom;
			else
				v(i) = vUpper;
			end
		end

	% CASE 2: no upper constraint, with lower constraint
	elseif ~isnan(vLower) && isnan(vUpper)
		for i = 1:length(v_in)
			if v_in(i) > vLower* denom
				v(i) = v_in(i) / denom;
			else
				v(i) = vLower;
			end
		end

	% CASE 3: no upper constraint, no lower constraint
	elseif isnan(vLower) && isnan(vUpper)
		for i = 1:length(v_in)
			v(i) = v_in(i) / denom;
		end

	% CASE 4: with upper constraint, with lower constraint
	elseif ~isnan(vLower) && ~isnan(vUpper)

		% check that constraints are appropriate
		if ~(vLower< vUpper)
			error('vLower must be less than vUpper')
		end

		lowerBound = vLower* denom;
		upperBound = vUpper* denom;

		% apply resolvent
		for i = 1:length(v_in)
			if v_in(i) > vLower* denom && v_in(i) < vUpper* denom
				v(i) = v_in(i) / denom;
			elseif v_in(i) <= vLower* denom
				v(i) = vLower;
			elseif v_in(i) >= vUpper* denom
				v(i) = vUpper;
			end
		end
	end

end
