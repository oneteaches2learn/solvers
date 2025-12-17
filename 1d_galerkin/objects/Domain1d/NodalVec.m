classdef NodalVec

    properties
        values  % vector of function values at mesh nodes
    end

    methods
        function self = NodalVec(vals)
        % NODALVEC Constructor for NodalVec object
        %
        % Inputs:
        %   dom      - Domain1d object
        %   vals     - function values at mesh nodes

            self.values = vals;

        end
    end
end