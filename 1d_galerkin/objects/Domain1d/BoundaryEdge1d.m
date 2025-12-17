classdef BoundaryEdge1d

    properties
		vertex
		outwardNormal
		boundaryType
		boundaryCondition
		boundaryCondition_ddu
		boundaryCondition_corr
		nodes
		nNodes
	end

    methods 
        % CONSTRUCTOR
        function self = BoundaryEdge1d(vertex)

            if nargin == 1
                self.vertex = vertex;
            end

        end


        % SETTERS
        function self = set.boundaryType(self,bcType)
            mustBeMember(bcType, {'D','N','R','P'});
            self.boundaryType = bcType;
        end

        function self = set.boundaryCondition(self,bc)
                   
            % if bc is cell, it represents a Robin condition and contains two functions
            if iscell(bc)
                for j = 1:length(bc)
                    bc{j} = Coefficients1d.func2functionHandle(bc{j}); % validate each function
                end
                self.boundaryCondition = bc;

            % else bc is a single function
            else
                bc = Coefficients1d.func2functionHandle(bc); % validate function
                self.boundaryCondition = bc;
            end
        end

        function self = set.boundaryCondition_ddu(self,bc_ddu)
                   
            % if bc_ddu is cell, it represents a Robin condition and contains two functions
            if iscell(bc_ddu)
                for j = 1:length(bc_ddu)
                    bc_ddu{j} = Coefficients1d.func2functionHandle(bc_ddu{j}); % validate each function
                end
                self.boundaryCondition_ddu = bc_ddu;

            % else bc_ddu is a single function
            else
                bc_ddu = Coefficients1d.func2functionHandle(bc_ddu); % validate function
                self.boundaryCondition_ddu = bc_ddu;
            end
        end


        % GETTERS
        function val = get.nNodes(self)
            val = length(self.nodes);   
        end

    end
end