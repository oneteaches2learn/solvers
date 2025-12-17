classdef AuxFun1d < Coefficients1d

    properties
        k     % diffusion coefficient, function handle
        r     % reaction coeficient, function handle
        dr_du % derivative of reaction coefficient wrt u, function handle
        f     % source term, function handle
    end

    methods
        % CONSTRUCTOR
        function self = AuxFun1d(k,f)

            % call superclass constructor
            self@Coefficients1d();

            % configure and storeinputs
            if nargin == 2

                self.k = Coefficients1d.func2functionHandle(k);
                self.f = Coefficients1d.func2functionHandle(f);

            end
        end


        % SETTERS 
        function self = set.k(self,k)
            self.k = Coefficients1d.func2functionHandle(k);
        end
        function self = set.r(self,r)
            self.r = Coefficients1d.func2functionHandle(r);
        end
        function self = set.dr_du(self,dr_du)
            self.dr_du = Coefficients1d.func2functionHandle(dr_du);
        end
        function self = set.f(self,f)
            self.f = Coefficients1d.func2functionHandle(f);
        end


        % GETTERS
        function cofs = cofs(self)
        % COFS returns a struct of the coefficients stored in the AuxFun1d object.
        %
        % NOTES:
        %  (1) This is meant to be called like a "getter" method. However, it is
        %  NOT a true getter, i.e. get.cofs, in that MATLAB does not support
        %  overwriting getters by subclasses. Since subclasses derived from
        %  AuxFun1d may want to add additional coefficients, we implement this
        %  as a regular method that can be called by subclasses and extended as
        %  needed.
        
            cofs.k = self.k;
            cofs.f = self.f;
            cofs.r = self.r;
            cofs.dr_du = self.dr_du;
        end

    end
end