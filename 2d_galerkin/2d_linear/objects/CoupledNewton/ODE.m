classdef ODE

    properties
        f
        u0
        t0
        T
        dt
        order
    end

    methods
        function obj = ODE()

            %{
            obj.f = f;
            obj.u0 = u0;
            obj.t0 = t0;
            obj.T = T;
            obj.dt = dt;
            obj.order = order;
            %}
        end

        function [t, u] = solve(obj)
            t = obj.t0:obj.dt:obj.T;
            u = zeros(length(t), length(obj.u0));
            u(1, :) = obj.u0;
            for i = 2:length(t)
                u(i, :) = obj.step(t(i), u(i-1, :));
            end
        end

        function u = step(obj, t, u)
            switch obj.order
                case 1
                    u = u + obj.dt * obj.f(t, u);
                case 2
                    u = u + obj.dt * obj.f(t, u) + 0.5 * obj.dt^2 * obj.f(t, u);
                otherwise
                    error('Invalid order');
            end
        end

    end
end