classdef KalmanFilter2d < CoupledNewtonSolver2d_rxndiff
% This class is made to hold special data and functions specific to the
% frostbite model. The "solve" method is called by the superclass constructor.

	properties
        ensemble
	end

	methods
		function self = KalmanFilter2d(dom,auxfun,ODE)

			% call superclass constructor
			self@CoupledNewtonSolver2d_rxndiff(dom,auxfun,ODE);

			% solve
			% ... called by superclass constructor

		end

    end
end
