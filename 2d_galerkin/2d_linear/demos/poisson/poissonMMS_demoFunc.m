function poissonMMS_demoFunc(params, options)
    % POISSONMMS_DEMO Function to run the Poisson MMS demo
    % Inputs:
    %   params  - Structure containing user input parameters
    %   options - Structure containing options including mmsparams

    % Extract user inputs from params structure
    xLim_dom = params.xLim_dom;
    yLim_dom = params.yLim_dom;
    xLim_Y = params.xLim_Y;
    yLim_Y = params.yLim_Y;
    incRatio = params.incRatio;
    eps = params.eps;
    bTypes_outer = params.bTypes_outer;
    bTypes_inner = params.bTypes_inner;
    k = params.k;
    r = params.r;
    demo = params.demo;
    uTrue = params.uTrue;

    % Extract mmsparams from options structure
    base = options.base;
    demo = options.demo;
    timeOffset = options.timeOffset;
    timeFactor = options.timeFactor;
    pmin = options.pmin;
    pmax = options.pmax;
    meshInclusions = options.meshInclusions;
    effectiveRegion = options.effectiveRegion;


    % MMS TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('MMS Test Begun\n')

    % Assemble inputs
    auxfun = ManufacturedFunctions2d_poisson(k, r, uTrue);
    mmsparams = MMSParams(base, ...
                            demo = options.demo, ...
                            timeOffset = options.timeOffset, ...
                            timeFactor = options.timeFactor, ...
                            pmin = options.pmin, ...
                            pmax = options.pmax, ...
				            meshInclusions = options.meshInclusions, ...
                            effectiveRegion = options.effectiveRegion);

    % Build dom_eps_epsain
    fprintf('Initialization\n')
    fprintf(' Constructing Domain:'), tic
        % dom = Domain2d(xLim_dom, yLim_dom);
        inc = Inclusion2d_circle(xLim_Y, yLim_Y, incRatio);
        % inc = Inclusion2d_square(xLim_Y, yLim_Y, incRatio);
        dom = Domain2d_punctured(xLim_dom, yLim_dom, inc, eps);
        dom = dom.add_yline(0.8);
        dom = dom.setBCTypes([bTypes_outer, bTypes_inner]);
        executionTime = toc;
    fprintf(' %f s\n', executionTime)

    % Run MMS test
    if demo == 0
        mms = GalerkinMMS2d_poisson(dom, auxfun, mmsparams, 'errType', 'L2');
    else
        mms = GalerkinMMS2d_poisson(dom, auxfun, mmsparams);
        prob = mms.problems{1};
    end
end