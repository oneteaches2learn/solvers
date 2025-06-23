classdef GalerkinFineGrid2d < GalerkinMMS2d

    properties
    end

    methods
        function self = GalerkinFineGrid2d(dom,auxfun,mmsparams,NameValueArgs)

            % call parent constructor
            self = self@GalerkinMMS2d();

			if nargin > 0

				% store inputs
                self.domain = dom;
				self.auxFunctions = auxfun;
				self.errType = "L2";
				self.mmsParams = mmsparams;
				try, time = NameValueArgs.time; end

				% run fine grid test
				self.problems = self.solveManufacturedProblems;

                % compute errors
                [self.errors,self.ratios,self.orders] = self.computeErrors;

			end

        end


		function problems = solveManufacturedProblems(self,varargin)

			% unpack variables
			cofs = self.auxFunctions;
			pmin = self.mmsParams.pmin;
			pmax = self.mmsParams.pmax;
			padd = self.mmsParams.pfine_add;

			% run mms test
			fprintf('Fine Grid Test Begun\n')
			fprintf('Solving Problems\n')
			ind = 1;
			for p = [pmin:pmax,pmax+padd]
				
				% successively refine mesh
				fprintf(' p = %i solved:',p); tic;
				dom_p = self.configureDomain(p); 
				prob_p = self.solve(dom_p,cofs);

				% store results
                problems{ind} = prob_p;

				% prepare next trial
				ind = ind + 1;
				executionTime = toc;
				fprintf(' %f s\n',executionTime)

			end
		end

		function dom_p = configureDomain(self,p)

			% unpack variables
			base = self.mmsParams.base;
			tOff = self.mmsParams.timeOffset;
			tFac = self.mmsParams.timeFactor;

            % generate path to mesh file
            root = self.mmsParams.mshRoot;
            filename = fullfile(root, sprintf('mesh%d.m', p));
            if ~isfile(filename)
                error('Mesh file does not exist: %s', filename);
            end

            % load mesh file; stores msh variable in workspace
            run(filename);

            % try to convert mesh to Domain2d object
            try 
                dom_p = Domain2d.domainFromGmsh(msh);

            % on error, try to reverse the orientation of triangles
            catch
                msh.TRIANGLES = msh.TRIANGLES(:,[1,3,2,4]);
                dom_p = Domain2d.domainFromGmsh(msh);
            end

            % update boundary
            dom_p.boundary = self.domain.boundary;
			dom_p = dom_p.setBoundaryNodes;

			% set time grid
			if isa(self,'GalerkinFineGrid2d_parabolic') || isa(self,'NewtonGalerkinMMS2d_parabolic')
				dom_p.time = self.domain.time;
				dom_p.time = dom_p.time.setMesh(tFac*(p-tOff),base);
			end
		end
        

		function [errors,ratios,orders] = computeErrors(self)

			% store variables
			trials = length(self.problems) - 1;
			base = self.problems{1}.domain.mesh.base;
			pmin = self.mmsParams.pmin;
			pmax = self.mmsParams.pmax;
			padd = self.mmsParams.pfine_add;
            quadType = self.mmsParams.quadType;

			% Store errors
			errors = zeros(1,trials);
			ratios = zeros(1,trials-1);
			orders = zeros(1,trials-1);

			% Compute error
			fprintf('Computing Errors:\n')

			for i = 1:trials

				fprintf(' Trial: '); tic;

				% if time varying, match time points
				if isa(self,'GalerkinFineGrid2d_parabolic')

					% find indices of fine grid times that match coarse grid times
					time_fine = find(ismember(self.problems{trials+1}.domain.time.tGrid, self.problems{i}.domain.time.tGrid));

				else 
					time_fine = 1;

				end

				% store fine solution on coarse solution time points
				uFine = self.problems{trials+1}.solution(:,time_fine);

				% NOTE: speedup can be attained if you vectorize how MidpointProjector works
				sol = zeros(size(uFine,1),length(time_fine));
				for j = 1:length(time_fine)

					% store coarse solution
					U_coarse = self.problems{i}.solution(:,j);

					% project coarse solution onto fine grid
					coarseLevel = pmin + (i - 1);
					fineLevel = pmax + padd;
					proj = MidpointProjector(U_coarse,coarseLevel,fineLevel);
					sol(:,j) = proj.U{fineLevel};
						
				end

				% compute quadrature (on each timestep, if time-varying)
				if strcmp(quadType,'threePoint')
					err = self.problems{trials+1}.domain.L2err_threePointQuadrature_nodal(sol,uFine);
				elseif strcmp(quadType,'centroid')
					err = self.problems{trials+1}.domain.L2err_centroidQuadrature_nodal(sol,uFine);
				elseif strcmp(quadType,'nodal')
					err = self.problems{trials+1}.domain.L2err_nodalQuadrature(sol,uFine);
				end

				% compute L_infty error across all time steps
				errors(i) = max(err);

				executionTime = toc;
				fprintf('%f s\n',executionTime) 
			end

			% Compute ratios and orders
			for i = 2:trials
				ratios(i-1) = errors(i-1)/errors(i);
				orders(i-1) = log(ratios(i-1)) / log(base);
			end

		end

    end

end



% GRAVEYARD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
		function dom_p = configureDomain(self,ind)

			% unpack variables
			%cofs = self.auxFunctions;
			base = self.mmsParams.base;
			pmin = self.mmsParams.pmin;
			pmax = self.mmsParams.pmax;
			tOff = self.mmsParams.timeOffset;
			tFac = self.mmsParams.timeFactor;
			region = self.mmsParams.effectiveRegion;
			inc_onoff = self.mmsParams.meshInclusions;
                
            % generate base mesh
            p = pmin;
            dom_p = self.domain;
            dom_p = dom_p.setMesh(p,base);

            % generate subsequent meshes by refining first mesh
            if ind ~= 1
                dom_p = dom_p.refineMesh(ind-1);
            end

			% assign boundary nodes to edges
			dom_p = dom_p.setBoundaryNodes;

			% set time grid
			if isa(self,'GalerkinMMS2d_parabolic') || isa(self,'NewtonGalerkinMMS2d_parabolic')
				dom_p.time = dom_p.time.setMesh(tFac*(p-tOff),base);
			end

		end
        %}
 %

            %{
			% check if mesh should be loaded from file
			try self.mmsParams.mshRoot ~= "";

                % generate domain from .m file
				root = self.mmsParams.mshRoot;
                filename = fullfile(root, sprintf('mesh%d.m', p));
                if ~isfile(filename)
                    error('Mesh file does not exist: %s', filename);
                end

                % load mesh file
                run(filename);
                try 
                    % try to load mesh
                    dom_p = Domain2d.domainFromGmsh(msh);
                catch
                    % on error, try to reverse the orientation of triangles
                    msh.TRIANGLES = msh.TRIANGLES(:,[1,3,2,4]);
                    dom_p = Domain2d.domainFromGmsh(msh);
                end

                % update boundary with stored boundary
                dom_p.boundary = self.domain.boundary;

			% otherwise generate mesh using MATLAB utilities
			catch 
				dom_p = self.domain;
				if isa(dom_p,'Domain2d_punctured')
					dom_p = dom_p.setMesh(p,base,meshInclusions=inc_onoff,effectiveRegion=region);
				else
					dom_p = dom_p.setMesh(p,base);
				end

			end
            %}

