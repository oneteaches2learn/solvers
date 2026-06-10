classdef StolwijkActivation
    % StolwijkActivation
    %
    % Function-handle factory for the Stolwijk--Hardy--Castellani
    % temperature-dependent activation coefficient.
    %
    % The primary method value() returns a function handle for
    %
    %     sigma(u,v) = (1 + c1*D(u,v)) / (1 + C(u,v)),
    %
    % where
    %
    %     D(u,v) = max(0, c2*(u-uStar) + c3*(v-vStar)),
    %     C(u,v) = max(0, -c4*(u-uStar) - c5*(v-vStar)).
    %
    % The methods du() and dv() return function handles for the
    % derivatives of sigma with respect to u and v.

    properties
        c1 double
        c2 double
        c3 double
        c4 double
        c5 double

        uStar double
        vStar double
    end

    methods
        function obj = StolwijkActivation(opts)
            arguments
                opts.c1 double = 0.1;
                opts.c2 double = 28424;
                opts.c3 double = 4870;
                opts.c4 double = 1.1;
                opts.c5 double = 3.3;
                opts.uStar double = 33;
                opts.vStar double = 37;
            end

            obj.c1 = opts.c1;
            obj.c2 = opts.c2;
            obj.c3 = opts.c3;
            obj.c4 = opts.c4;
            obj.c5 = opts.c5;
            obj.uStar = opts.uStar;
            obj.vStar = opts.vStar;
        end


        % FUNCTION HANDLES
        function fh = value(obj)
            % VALUE Return function handle for sigma(u,v).

            c1 = obj.c1;
            D = obj.D;
            C = obj.C;

            fh = @(u,v)( (1 + c1 .* D(u,v)) ./ (1 + C(u,v)) );

        end

        function fh = D(obj)
            % D Return function handle for the vasodilatory term D(u,v).

            c2 = obj.c2;
            c3 = obj.c3;
            uStar = obj.uStar;
            vStar = obj.vStar;

            fh = @(u,v) 1e-6 * max(0, c2 .* (v - vStar) + c3 .* (u - uStar));
        end

        function fh = C(obj)
            % C Return function handle for the vasoconstrictive term C(u,v).

            c4 = obj.c4;
            c5 = obj.c5;
            uStar = obj.uStar;
            vStar = obj.vStar;

            fh = @(u,v) max(0, -c4 .* (v - vStar) - c5 .* (u - uStar));
        end


        % DERIVATIVE FUNCTIONS
        function fh = du(obj)
            % DU Return function handle for d/du sigma(u,v).
            %
            % At the max() kink, the derivative is chosen to be zero.

            c1 = obj.c1;
            D    = obj.D();
            C    = obj.C();
            D_du = obj.D_du();
            C_du = obj.C_du();

            fh = @(u,v) ...
                (c1 .* D_du(u,v) .* (1 + C(u,v)) ...
                - (1 + c1 .* D(u,v)) .* C_du(u,v)) ...
                ./ (1 + C(u,v)).^2;
        end

        function fh = D_du(obj)
            % D_du returns the u-partial derivative of D(u,v)

            c2 = obj.c2;
            c3 = obj.c3;
            uStar = obj.uStar;
            vStar = obj.vStar;

            fh = @(u,v) 1e-6 * c2 .* double(c2 .* (v - vStar) + c3 .* (u - uStar) > 0);

        end
        
        function fh = C_du(obj)
            % C_du returns the u-partial derivative of C(u,v)

            c4 = obj.c4;
            c5 = obj.c5;
            uStar = obj.uStar;
            vStar = obj.vStar;

            fh = @(u,v) -c4 .* double(-c4 .* (u - uStar) - c5 .* (v - vStar) > 0);

        end
    end
end