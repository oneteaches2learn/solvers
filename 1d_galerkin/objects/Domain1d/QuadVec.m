classdef QuadVec

    properties
        values  % vector of function values at quadrature points
    end

    methods
        function self = QuadVec(vals)
        % QUADVEC Constructor for QuadVec object
        %
        % Inputs:
        %   vals     - function values at quadrature points

            self.values = vals;

        end
    end
end