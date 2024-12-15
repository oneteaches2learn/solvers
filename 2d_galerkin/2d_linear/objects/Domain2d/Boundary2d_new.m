classdef Boundary2d_new

    properties
        effectiveRegion = "Omega_eps"
        dl
    end

    properties (Dependent)
        nEdges
        edges
    end

    properties (Dependent,Hidden)
        nEdges_outer
        nEdges_inclusions
    end

    properties (Hidden)
    end


    methods
        function self = Boundary2d_new(dl)

            self.dl = dlObj(dl);

        end

        % GETTERS: number of edges
        function val = get.nEdges_outer(self)

            val = 4;

        end

        function val = get.nEdges_inclusions(self)

            if ~strcmp(self.effectiveRegion,"Omega_eps")
                val = 0;
            else
                val = self.dl.nSeg_inclusions;
            end

        end

        function val = get.nEdges(self)

            val = self.nEdges_outer + self.nEdges_inclusions;

        end

		function self = update_outerSegmentIDs(self)

			% update outer edge IDs
			%n_outerSegments = size(self.dl_outer,2)
			%outerIDs = 1:n_outerSegments;
			%outerIDs = reshape(outerIDs,[],2)';

			%southIDs = outerIDs(1,1);
			%eastIDs  = outerIDs(1,2:end);
			%northIDs = outerIDs(2,1);
			%westIDs  = outerIDs(2,2:end);

			%self.edges(1).segmentIDs = southIDs;
			%self.edges(2).segmentIDs = eastIDs;
			%self.edges(3).segmentIDs = northIDs;
			%self.edges(4).segmentIDs = westIDs;

		end

		function self = update_inclusionSegmentIDs(self)

			% get number of outer segments
			%n_outerSegments = size(self.dl_outer,2);

			% update inclusion IDs
			%if length(self.edges) > 4
			%	for i = 5:length(self.edges)
			%		j = i - 4;
			%		self.edges(i).segmentIDs = j + n_outerSegments;
			%	end
			%end

		end

        %{
		function self = add_yline(self,varargin)

			% store variables
			dl_outer = self.dl.dl_outer;
			dl_inclusions = self.dl.dl_inclusions;

			% collect values where yLines should be added
			if nargin == 1, yBar = mean(self.dl.yLim);
			else, for i = 1:length(varargin), yBar(i) = varargin{i}; end, end
			
			% build new dl matrix
			dl_outer = self.split_dl_outer(yBar,dl_outer);
			dl_outer = self.increment_dl_outer(yBar,dl_outer);
			dl_inclusions = self.increment_dl_inclusions(yBar,dl_inclusions);
			dl_yLines = self.build_dl_yLines(dl_outer);
			
			% store dl matrix
			self.dl = dlObj([dl_outer,dl_inclusions,dl_yLines]);

			% update edge IDs
			self = self.update_dl_IDs(dl_outer,dl_inclusions,dl_yLines);
			self = self.update_outerSegmentIDs;
			self = self.update_inclusionSegmentIDs;

			pdegplot([dl_outer, dl_inclusions, dl_yLines],'EdgeLabels','on','FaceLabels','on')

		end
        %}

        function self = add_yline(self,varargin)

            class(self.dl_new)
            % update dl
            %self.dl = self.dl_new.add_yline(varargin{:});



            %pdegplot(self.dl.mat,'EdgeLabels','on','FaceLabels','on')

        end

    end
end