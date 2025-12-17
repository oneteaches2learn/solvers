classdef AuxFun1d_projection < AuxFun1d

    properties
        lambda    % penalty parameter
        u_fixed   % known function u, function handle
        v_fixed   % known function v, function handle
    end

    properties (Hidden)
    end

    methods
        % CONSTRUCTOR
        function self = AuxFun1d_projection(k,f)

            % call superclass constructor
            self@AuxFun1d(k,f);

        end

        % SETTERS
        function self = set.lambda(self,lambda)
            self.lambda = lambda;
        end 

        function self = set.u_fixed(self,u_fixed)
            self.u_fixed = self.func2functionHandle(u_fixed);
        end

        function self = set.v_fixed(self,v_fixed)
            self.v_fixed = self.func2functionHandle(v_fixed);
        end

        % GETTERS
        function cofs = cofs(self)
            
            % call superclass method
            cofs = cofs@AuxFun1d(self);

            % add on coefficients
            cofs.lambda = self.lambda;
            cofs.u_fixed = self.u_fixed;
            cofs.v_fixed = self.v_fixed;
        end

    end
end