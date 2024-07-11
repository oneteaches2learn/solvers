classdef Transmissibility

	properties
		cellEdges
	end

	methods
		function obj = Transmissibility(domain,permeability)
			obj.check_inputs(domain,permeability);
			obj.cellEdges = obj.compute_transmissibility(domain,permeability);
		end

		function check_inputs(obj,domain,permeability)
			err = "Transmissibility inputs must be a Domain1d object and a Permeablity object";
			if ~isa(domain,'Domain') || ~isa(permeability,'Permeability')
				error(err)
			end
		end

		function tx = compute_transmissibility(obj,domain,permeability)
		%COMPUTETRANSMISSIBILITY computes the transmissibility at each edge

			dx = domain.cellWidths;
		    n  = domain.cellNum;
			perm = permeability.cellEdges;

		    tx = zeros(1,n+1);

			tx(2:n) = obj.harmonic_mean(perm(2:n)./dx(1:n-1),perm(2:n)./dx(2:n));
			tx(1) = 2*perm(1)/dx(1);
		    tx(n+1) = 2*perm(n+1)/dx(n);
		end

		function val = harmonic_mean(obj,x,y)
		%HARMONICMEAN(x,y) computes harmonic mean of x and y
    
		    val = 2./(1./x + 1./y);

		end

	end
end


