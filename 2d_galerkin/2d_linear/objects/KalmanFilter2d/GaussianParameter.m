classdef GaussianParameter
    properties
        mean double
        stdDev double
    end

    methods
        function obj = GaussianParameter(meanValue, stdDev)

            if any(stdDev(:) < 0)
                error("GaussianParameter:InvalidStdDev", ...
                      "Standard deviation must be nonnegative.");
            end

            if ~isequal(size(meanValue), size(stdDev))
                error("GaussianParameter:SizeMismatch", ...
                    "meanValue and stdDev must have the same size.");
            end

            obj.mean = meanValue;
            obj.stdDev = stdDev;
        end

        function x = sample(obj, varargin)
            if isempty(varargin)
                sampleSize = size(obj.mean);
            else
                sampleSize = cell2mat(varargin);
            end

            x = obj.mean + obj.stdDev .* randn(sampleSize);
        end

        function objSlice = subsref(obj, S)
            if strcmp(S(1).type, "()")
                objSlice = GaussianParameter( ...
                    obj.mean(S(1).subs{:}), ...
                    obj.stdDev(S(1).subs{:}));

                if numel(S) > 1
                    objSlice = builtin("subsref", objSlice, S(2:end));
                end
            else
                objSlice = builtin("subsref", obj, S);
            end
        end

        function x = sampleTruncated(obj, lower, upper)
            maxTries = 10000;

            for k = 1:maxTries
                x = obj.sample();

                if all(x(:) >= lower) && all(x(:) <= upper)
                    return
                end
            end

            error("GaussianParameter:TruncationFailed", ...
                "Failed to sample inside bounds after many attempts.");
        end

        function x = samplePositive(obj)
            x = obj.sampleTruncated(0, Inf);
        end

        function v = variance(obj)
            v = obj.stdDev.^2;
        end
    end

    methods (Static)
        function obj = fromVariance(meanValue, variance)
            if any(variance(:)) < 0
                error("GaussianParameter:InvalidVariance", ...
                      "Variance must be nonnegative.");
            end

            obj = GaussianParameter(meanValue, sqrt(variance));
        end
    end
end