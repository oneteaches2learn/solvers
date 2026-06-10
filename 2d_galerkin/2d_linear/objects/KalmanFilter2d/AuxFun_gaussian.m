classdef AuxFun_gaussian
    % AuxFun_gaussian
    %
    % Stores frostbite-style auxiliary functions with possibly-random
    % Gaussian parameters. Calling sample() produces a deterministic auxfun
    % suitable for passing to an existing solver.
    %
    % Current frostbite-style inputs:
    %   c      heat capacity coefficient
    %   k      conductivity coefficient
    %   f      source term
    %   uInit  initial condition
    %
    % For now, r is set to zero.

    properties
        c
        k
        f
        uInit
    end

    methods
        function obj = AuxFun_gaussian(c, k, f, uInit)
            obj.c = c;
            obj.k = k;
            obj.f = f;
            obj.uInit = uInit;
        end

        function [auxfun, samples] = sample(obj)
            % SAMPLE
            %
            % Produces a deterministic auxfun realization.
            %
            % Usage:
            %   auxfun = obj.sample();
            %   [auxfun, samples] = obj.sample();

            samples = struct();

            [cSample, samples.c]         = obj.sampleEntry(obj.c);
            [kSample, samples.k]         = obj.sampleEntry(obj.k);
            [fSample, samples.f]         = obj.sampleEntry(obj.f);
            [uInitSample, samples.uInit] = obj.sampleEntry(obj.uInit);

            auxfun = struct();
            auxfun.cofs = struct();

            auxfun.cofs.c = obj.toSpatialFunctionHandle(cSample);
            auxfun.cofs.k = obj.toSpatialFunctionHandle(kSample);
            auxfun.f      = obj.toSpatialFunctionHandle(fSample);
            auxfun.uInit  = obj.toSpatialFunctionHandle(uInitSample);

            % Frostbite placeholder. We are setting aside r for now.
            auxfun.cofs.r = @(x1, x2, t, u, v) zeros(size(u));
            auxfun.cofs.dr_du = @(x1, x2, t, u, v) zeros(size(u));
        end

        function [auxfuns, samples] = sampleEnsemble(obj, N)
            % SAMPLEENSEMBLE
            %
            % Produces N deterministic auxfun realizations.
            %
            % Usage:
            %   auxfuns = obj.sampleEnsemble(N);
            %   [auxfuns, samples] = obj.sampleEnsemble(N);

            auxfuns = cell(N, 1);
            samples = repmat(struct(), N, 1);

            for j = 1:N
                [auxfuns{j}, samples(j)] = obj.sample();
            end
        end
    end

    methods (Access = private)
        function [valueSampled, sampleRecord] = sampleEntry(obj, value)
            % SAMPLEENTRY
            %
            % If value has a sample() method, call it.
            % Otherwise, leave value unchanged.

            %#ok<INUSD>

            if obj.hasSampleMethod(value)
                valueSampled = value.sample();
                sampleRecord = valueSampled;
            else
                valueSampled = value;
                sampleRecord = value;
            end
        end

        function tf = hasSampleMethod(obj, value)
            %#ok<INUSD>

            tf = false;

            if isobject(value)
                try
                    tf = ismethod(value, "sample");
                catch
                    tf = false;
                end
            end
        end

        function fh = toSpatialFunctionHandle(obj, value)
            % TOSPATIALFUNCTIONHANDLE
            %
            % Converts numeric, symbolic, symfun, or function_handle inputs
            % into a spatial function handle of the form
            %
            %   fh(x1, x2)
            %
            % This matches the current parabolic/frostbite coefficient style.

            %#ok<INUSD>

            if isa(value, "function_handle")
                fh = value;
                return
            end

            if isnumeric(value)
                fh = @(x1, x2) value + zeros(size(x1));
                return
            end

            if isa(value, "symfun")
                fh = matlabFunction(value);
                return
            end

            if isa(value, "sym")
                x = sym("x", [1 2], "real");
                value = symfun(value, x);
                fh = matlabFunction(value);
                return
            end

            error("AuxFun_gaussian:InvalidInput", ...
                  "Input must be numeric, symbolic, symfun, function_handle, or an object with sample().");
        end
    end
end