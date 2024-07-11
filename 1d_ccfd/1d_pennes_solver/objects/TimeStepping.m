classdef TimeStepping

	properties
		stepSize
		stepNum
		start
		stop
	end

	methods
		function obj = TimeStepping(stepSize,stepNum,start)
			obj.stepSize = stepSize;
			obj.stepNum = stepNum;

			if nargin == 3
				obj.start = start;
			else
				obj.start = 0;
			end

			obj.stop = obj.start + stepSize * stepNum;
		end
	end
end

