classdef MidpointProjector
% MIDPOINTPROJECTOR - Project a solution from a coarse mesh to a fine mesh
%   This class is used to project a solution from a coarse mesh to a fine
%   mesh. The projection is done by first finding the shared nodes between
%   the two meshes and assigning the solution values at those nodes. Then,
%   the solution is projected onto the midpoints of the edges in the
%   coarse mesh. The average of the solution values at the endpoints of
%   the edges is assigned to the midpoint nodes in the fine mesh.
%

    properties
        startingLevel   % int, refinement level of the coarsest mesh
        endingLevel     % int, refinement level of the finest mesh
        U               % cell, solution vectors 
    end

    properties (Hidden)
        U_coarse        % double, Solution on the coarse mesh
        U_fine          % double, Projected solution on the fine mesh
    end

    methods
        function self = MidpointProjector(U_coarse,startingLevel,endingLevel)
        % Constructor for the MidpointProjector class

            if nargin == 3
                % Store inputs
                self.startingLevel = startingLevel;
                self.endingLevel = endingLevel;
                self.U{startingLevel,1} = U_coarse;
                self.U_coarse = U_coarse;

                self = self.project;

            end


        end

        function self = project(self)

            % loop over starting levels
            for i = self.startingLevel:self.endingLevel-1

                % project onto the next mesh
                self = self.project_step(i);

                % update system state
                self.U_coarse = self.U{i+1,1};

            end
            
        end

        function self = project_step(self,currentLevel)
            
            % store domains
            dom_coarse = self.getDom(currentLevel);
            dom_fine = self.getDom(currentLevel+1);

            % instantiate projection vector
            U_coarse = self.U_coarse;
            U_fine = zeros(dom_fine.mesh.nNodes,1);

            % locate coarse nodes in fine mesh
            [found, idx_shared] = ismember( ...
                single(dom_coarse.mesh.nodes), ...
                single(dom_fine.mesh.nodes), ...
                'rows');
            if ~all(found)
                error('Some coarse nodes not found in fine mesh nodes.');
            end

            % store the solution values at the shared nodes
            U_fine(idx_shared) = U_coarse; 

            % find indices of all coarse midpoints in fine mesh nodes
            [found, idx_midpts] = ismember( ...
                single(dom_coarse.mesh.midpoints), ...
                single(dom_fine.mesh.nodes), ...
                'rows' );
            if ~all(found)
                error('Some midpoints not found in the fine mesh nodes.');
            end    

            % interpolate the coarse solution on the midpoints of the edges
            sol_coarse_midpts = (U_coarse(dom_coarse.mesh.edges(:,1),:) + ...
                                    U_coarse(dom_coarse.mesh.edges(:,2),:)) / 2;

            % assign the average solution value to appropriate location in U_1proj2
            U_fine(idx_midpts) = sol_coarse_midpts;

            % store the projected solution
            self.U{currentLevel+1,1} = U_fine;

        end

        function dom = getDom(self,requested_level)

            % load mesh
            dir = '/home/tyler/Software/MATLAB/figures/2025_amc_spring_talk/coupled_mms/meshes';
            filename = sprintf('regular_mesh%d.m',requested_level);
            fid = fullfile(dir,filename);
            run(fid);

            % reverse orientation of triangles
            msh.TRIANGLES = msh.TRIANGLES(:,[1 3 2]);

            % construct domain
            dom = Domain2d.domainFromGmsh(msh);

        end


        function plot(self,level)
        % plot plots the mesh at the given refinement level
        %   plot takes advantage of the existing plot function of the
        %   GalerkinSolver2d class.
        
            prob = GalerkinSolver2d();
            prob.domain = self.getDom(level);
            prob.solution = self.U{level};
            prob.plot(1);

        end

    end
end
 
