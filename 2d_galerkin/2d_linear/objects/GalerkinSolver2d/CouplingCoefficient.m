classdef CouplingCoefficient
    % CouplingCoefficient
    %
    % Represents the coupled PDE--ODE exchange terms
    %
    %     r(u,v) = A_b * sigma(u,v) * (u - v),
    %
    %     s(<u>,v) = chiVol * A_b * sigma(<u>,v) * (v - <u>),
    %
    % where u is the hand temperature, v is the core temperature, and
    % sigma is an activation law such as StolwijkActivation. The activation
    % coefficient sigma is dimensionless and models temperature-dependent
    % vasoconstriction/vasodilation.
    %
    % The coefficient A_b is the baseline blood--tissue exchange coefficient
    %
    %     A_b = c_b * rho_b * omega_b.
    %
    % Using the default values from the thesis,
    %
    %     c_b     = 3617      J kg^{-1} degC^{-1},
    %     rho_b   = 1050      kg m^{-3},
    %     omega_b = 1.1e-3    s^{-1},
    %
    % gives
    %
    %     A_b = 4177.635      W m^{-3} degC^{-1}.
    %
    % The factor chiVol accounts for the relative size of the hand and core
    % compartments in the ODE exchange term.

    properties
        Ab double
        chiVol double
        sigma
    end

    methods
        function obj = CouplingCoefficient(opts)

            arguments
                opts.Ab double = 4177.635;
                opts.chiVol double = 0.02;
                opts.sigma = StolwijkActivation();
            end

            obj.Ab = opts.Ab;
            obj.chiVol = opts.chiVol;
            obj.sigma = opts.sigma;
        end

        function y = r(obj, u, v)

            sigma = obj.activation.valueHandle();
            AbLocal = obj.Ab;

            fh = @(x1,x2,t,u,v) ...
                AbLocal .* sigma(u,v) .* (u - v);

            %y = obj.Ab .* obj.sigma.value(u,v) .* (u - v);
        end

        function y = dr_du(obj, u, v)
            sig = obj.sigma.value(u,v);
            sig_u = obj.sigma.du(u,v);

            y = obj.Ab .* (sig_u .* (u - v) + sig);
        end

        function y = dr_dv(obj, u, v)
            sig = obj.sigma.value(u,v);
            sig_v = obj.sigma.dv(u,v);

            y = obj.Ab .* (sig_v .* (u - v) - sig);
        end

        function y = s(obj, uAvg, v)
            y = obj.chiVol .* obj.Ab .* obj.sigma.value(uAvg,v) .* (v - uAvg);
        end

        function y = ds_du(obj, uAvg, v)
            sig = obj.sigma.value(uAvg,v);
            sig_u = obj.sigma.du(uAvg,v);

            y = obj.chiVol .* obj.Ab .* ...
                (sig_u .* (v - uAvg) - sig);
        end

        function y = ds_dv(obj, uAvg, v)
            sig = obj.sigma.value(uAvg,v);
            sig_v = obj.sigma.dv(uAvg,v);

            y = obj.chiVol .* obj.Ab .* ...
                (sig_v .* (v - uAvg) + sig);
        end
    end
end