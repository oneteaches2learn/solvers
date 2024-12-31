classdef ODE

    properties
        f
        vInit
        s
        t0
        T
        dt
        order
        solution
    end

    methods
        % CONSTRUCTOR
        function self = ODE(data,options)

            self.f = data.f;
            self.vInit = data.vInit;
            self.s = data.cofs.s;
            self.T = data.time.T;
            self.dt = data.time.dt;
            self.order = options.order;

        end

        function [t, u] = solve(obj)

            t = self.t0:self.dt:self.T;
            u = zeros(length(t), length(self.u0));
            u(1, :) = self.u0;
            for i = 2:length(t)
                u(i, :) = self.step(t(i), u(i-1, :));
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