classdef NodeObservationOperator
    properties
        obsNodes (:,1) double
        stateDimension (1,1) double
    end

    properties (Dependent)
        numObservations
        matrix
    end

    methods
        function obj = NodeObservationOperator(obsNodes, stateDimension)

            arguments
                obsNodes (:,1) double = [433; 2097; 5077]
                stateDimension (1,1) double = 5077
            end

            obj.obsNodes = obsNodes(:);
            obj.stateDimension = stateDimension;
        end

        function m = get.numObservations(obj)
            m = numel(obj.obsNodes);
        end

        function M = get.matrix(obj)
            m = obj.numObservations;
            M = sparse(1:m, obj.obsNodes, ones(m,1), ...
                       m, obj.stateDimension);
        end

        function y = observe(obj, psi)
            y = psi(obj.obsNodes);
        end
    end
end