function [pass,mms] = pennesMMS_test(dom,params,options)

	% MMS TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	fprintf('MMS Test Begun\n')

	% assemble inputs
	auxfun    = ManufacturedFunctions2d_pennes(params.c,params.k,params.r,params.uStar,params.uTrue);
	mmsparams = MMSParams(options.base, ...
							demo = options.demo, ...
							timeOffset = options.timeOffset, ...
							timeFactor = options.timeFactor, ...
							pmin = options.pmin, ...
							pmax = options.pmax, ...
							meshInclusions = options.meshInclusions, ...
							effectiveRegion = options.effectiveRegion);

	% run mms test
	mms = GalerkinMMS2d_pennes(dom,auxfun,mmsparams,errType="L2");


	% ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Check that the rate is converging toward the expected rate
	absErr = abs(mms.orders - 2);
	cond1 = all(absErr(1:end-1) > absErr(2:end));
	
	% Check that last reported rate is within tolerance of expected rate
	cond2 = abs(mms.orders(end) - options.expectedRate) < options.convergenceTolerance;

	% Check if error is extremely small
	cond3 = mms.errors(end) < 1e-12;
	
	% Use cond2 to determine pass/fail
	if cond2 || cond3
		pass = 1;
	else
		pass = 0;
	end

end