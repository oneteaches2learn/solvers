classdef EnsembleState
    % EnsembleState
    %
    % Stores an ensemble matrix and computes basic EnKF statistics.
    %
    % states(:,j) is the j-th ensemble member.
    % If states is nState-by-N, then
    %
    %   mean        = average over ensemble members
    %   perturbations = states - mean
    %   covariance = perturbations * perturbations' / (N - 1)

    properties
        states double   % nState-by-N ensemble matrix
    end

    properties (Dependent)
        numMembers      % Number of ensemble members
        stateDimension  % Dimension of each state vector
        ensembleMean    % nState-by-1 ensemble mean
        perturbations   % nState-by-N centered ensemble matrix
        covariance      % nState-by-nState sample covariance matrix
        variance        % Componentwise sample variance
        stdDev          % Componentwise sample standard deviation
    end

    methods
        function obj = EnsembleState(states)
            arguments
                states double
            end

            obj.states = states;
        end

        function obj = set.states(obj, states)
            if ndims(states) ~= 2
                error("EnsembleState:InvalidStates", ...
                      "states must be a 2D matrix.");
            end

            if size(states, 2) < 2
                error("EnsembleState:NotEnoughMembers", ...
                      "The ensemble must contain at least two members.");
            end

            obj.states = states;
        end

        function N = get.numMembers(obj)
            N = size(obj.states, 2);
        end

        function n = get.stateDimension(obj)
            n = size(obj.states, 1);
        end

        function m = get.ensembleMean(obj)
            m = mean(obj.states, 2);
        end

        function X = get.perturbations(obj)
            X = obj.states - obj.ensembleMean;
        end

        function C = get.covariance(obj)
            X = obj.perturbations;
            C = (X * X.') / (obj.numMembers - 1);
        end

        function v = get.variance(obj)
            X = obj.perturbations;
            v = sum(X.^2, 2) / (obj.numMembers - 1);
        end

        function s = get.stdDev(obj)
            s = sqrt(obj.variance);
        end
    end
end