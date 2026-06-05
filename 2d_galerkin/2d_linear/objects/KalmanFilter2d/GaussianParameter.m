classdef GaussianParameter
    properties
        mean double
        stdDev double
    end

    methods
        function obj = GaussianParameter(meanValue, stdDev)
            if stdDev < 0
                error("GaussianParameter:InvalidStdDev", ...
                      "Standard deviation must be nonnegative.");
            end

            obj.mean = meanValue;
            obj.stdDev = stdDev;
        end

        function x = sample(obj, varargin)
            if isempty(varargin)
                sampleSize = [1, 1];
            else
                sampleSize = cell2mat(varargin);
            end

            x = obj.mean + obj.stdDev * randn(sampleSize);
        end

        function x = sampleTruncated(obj, lower, upper)
            maxTries = 10000;

            for k = 1:maxTries
                x = obj.sample();
                if x >= lower && x <= upper
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
            v = obj.stdDev^2;
        end
    end

    methods (Static)
        function obj = fromVariance(meanValue, variance)
            if variance < 0
                error("GaussianParameter:InvalidVariance", ...
                      "Variance must be nonnegative.");
            end

            obj = GaussianParameter(meanValue, sqrt(variance));
        end
    end
end