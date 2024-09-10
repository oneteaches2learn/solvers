classdef Boundary2d_punctured < Boundary2d

    
    methods
        % CONSTUCTOR
        function self = Boundary2d_punctured(dl,bcTypes,bcConds)
            
            if nargin < 2
                bcTypes = 'skip';
            end
            if nargin < 3
                bcConds = 'skip';
            end

            % call superclass constructor
            self = self@Boundary2d(dl,bcTypes,bcConds);
            
        end

		function nEdges = set_nEdges(self,effectiveRegion)

            if strcmp(effectiveRegion,'Omega_eps')
                nEdges = size(self.dl_outer,2) + size(self.dl_inclusions,2);
            else
                nEdges = size(self.dl_outer,2);
            end

		end
    end
end