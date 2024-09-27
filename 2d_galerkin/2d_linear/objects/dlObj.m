classdef dlObj

    properties
        mat
    end

    properties (Dependent)
        nSeg_total
        nSeg_outer
        nSeg_inclusions
        nSeg_yLines
        segIDs_outer
        segIDs_inclusions
        segIDs_yLines
        dl_outer
        dl_inclusions
        dl_yLines
        segEdgeDict
        edgeSegDict
        edgeSegDict_outer
        edgeSegDict_inclusions
        edgeCoords
        edgeCoords_outer
        edgeCoords_inclusions
        faceIDs_all
        faceIDs_Omega_eps
        faceIDs_Q_eps
    end

    properties (Dependent,Hidden)
        xLowVal
        xHighVal
        yLowVal
        yHighVal
        xLim
        yLim
        edgeSegDict_outerSouth
        edgeSegDict_outerEast
        edgeSegDict_outerNorth
        edgeSegDict_outerWest
        edgeCoords_outerSouth
        edgeCoords_outerEast
        edgeCoords_outerNorth
        edgeCoords_outerWest
    end

    properties (Hidden,Constant)
        typeRow = 1
        xLowRow = 2
        xHighRow = 3
        yLowRow = 4
        yHighRow = 5
        leftRow = 6
        rightRow = 7
        centerxRow = 8
        centeryRow = 9
        radiusRow = 10
    end

    methods
        % CONSTRUCTOR
        function self = dlObj(mat)
            self.mat = self.set_dl(mat);
        end

        % SETTERS
        function mat = set_dl(self,mat)

            nRows = size(mat,1);
            if nRows ~= 10
                nCols = size(mat,2);
                pad = zeros(10 - nRows,nCols);
                mat = [mat; pad];
            end
        end

        % GETTERS: number of segments
        function val = get.nSeg_total(self)

            val = size(self.mat,2);

        end

        function val = get.nSeg_outer(self)
        
            val = sum(self.mat(self.rightRow,:) == 0);
         
        end

        function val = get.nSeg_yLines(self)

            val = sum( ...
            self.mat(self.xLowRow,:) == self.xLowVal & ...
            self.mat(self.xHighRow,:) == self.xHighVal & ...
            self.mat(self.yLowRow,:) > self.yLowVal & ...
            self.mat(self.yHighRow,:) < self.yHighVal);

        end

        function val = get.nSeg_inclusions(self)

            val = self.nSeg_total - self.nSeg_outer - self.nSeg_yLines;

        end

        % GETTERS: domain boundaries
        function val = get.xLowVal(self)
            val = min(self.mat(self.xLowRow:self.xHighRow,:),[],'all');
        end

        function val = get.xHighVal(self)
            val = max(self.mat(self.xLowRow:self.xHighRow,:),[],'all');
        end

        function val = get.yLowVal(self)
            val = min(self.mat(self.yLowRow:self.yHighRow,:),[],'all');
        end

        function val = get.yHighVal(self)
            val = max(self.mat(self.yLowRow:self.yHighRow,:),[],'all');
        end

        function val = get.xLim(self)
            val = [self.xLowVal self.xHighVal];
        end

        function val = get.yLim(self)
            val = [self.yLowVal self.yHighVal];
        end

        % GETTERS: segment IDs
        function val = get.segIDs_outer(self)
        
            val = find(self.mat(self.rightRow,:) == 0);
         
        end

        function val = get.segIDs_yLines(self)

            val = find( ...
            self.mat(self.xLowRow,:) == self.xLowVal & ...
            self.mat(self.xHighRow,:) == self.xHighVal & ...
            self.mat(self.yLowRow,:) > self.yLowVal & ...
            self.mat(self.yHighRow,:) < self.yHighVal);

        end

        function val = get.segIDs_inclusions(self)

            segIDs_total = 1:self.nSeg_total;
            temp = setdiff(segIDs_total,self.segIDs_outer);
            val = setdiff(temp,self.segIDs_yLines);

        end

        % GETTERS: matrix partitions
        function val = get.dl_outer(self)

            val = self.mat(:,self.segIDs_outer);

        end

        function val = get.dl_inclusions(self)

            val = self.mat(:,self.segIDs_inclusions);

        end

        function val = get.dl_yLines(self)

            val = self.mat(:,self.segIDs_yLines);

        end


        % GETTERS: segment-edge dictionary
        function val = get.segEdgeDict(self)

            val = [];
            for i = 1:length(self.edgeSegDict)
                
                val = [val i * ones(1,length(self.edgeSegDict{i}))];
                    
            end

        end


        % GETTERS: edge-segment dictionaries
        function val = get.edgeSegDict_outer(self)

            val{1} = self.edgeSegDict_outerSouth;
            val{2} = self.edgeSegDict_outerEast;
            val{3} = self.edgeSegDict_outerNorth;
            val{4} = self.edgeSegDict_outerWest;

        end

        function val = get.edgeSegDict_inclusions(self)

            if self.nSeg_inclusions == 0, val = {};
            else, for i = 1:length(self.segIDs_inclusions)
                val{i} = i + self.nSeg_outer;
            end, end

        end

        function val = get.edgeSegDict(self)
            
            val = [self.edgeSegDict_outer self.edgeSegDict_inclusions];

        end

        function val = get.edgeSegDict_outerSouth(self)

            val = find( ...
            self.dl_outer(self.yLowRow,:) == self.yLowVal & ...
            self.dl_outer(self.xLowRow,:) ~= self.dl_outer(self.xHighRow,:));

        end

        function val = get.edgeSegDict_outerEast(self)

            val = find( ...
            self.dl_outer(self.xHighRow,:) == self.xHighVal & ...
            self.dl_outer(self.yLowRow,:) ~= self.dl_outer(self.yHighRow,:));

        end

        function val = get.edgeSegDict_outerNorth(self)

            val = find( ...
            self.dl_outer(self.yHighRow,:) == self.yHighVal & ...
            self.dl_outer(self.xLowRow,:) ~= self.dl_outer(self.xHighRow,:));

        end

        function val = get.edgeSegDict_outerWest(self)

            val = find( ...
            self.dl_outer(self.xLowRow,:) == self.xLowVal & ...
            self.dl_outer(self.yLowRow,:) ~= self.dl_outer(self.yHighRow,:));

        end

        % GETTERS: face IDs
        function val = get.faceIDs_all(self)

            val = [self.faceIDs_Omega_eps self.faceIDs_Q_eps];

        end

        function val = get.faceIDs_Omega_eps(self)

            val = unique(self.dl_outer(self.leftRow,:));

        end

        function val = get.faceIDs_Q_eps(self)

            val = unique(self.dl_inclusions(self.leftRow,:));
            if val == 0, val = []; end

        end

        % GETTERS: edge coordinates
        function val = get.edgeCoords_outerSouth(self)

            val = [self.dl_outer(self.xLowRow,self.edgeSegDict_outerSouth(1)), ...
                   self.dl_outer(self.yLowRow,self.edgeSegDict_outerSouth(1));
                   self.dl_outer(self.xHighRow,self.edgeSegDict_outerSouth(end)), ...
                   self.dl_outer(self.yHighRow,self.edgeSegDict_outerSouth(end))];

        end

        function val = get.edgeCoords_outerEast(self)

            val = [self.dl_outer(self.xLowRow,self.edgeSegDict_outerEast(1)), ...
                   self.dl_outer(self.yLowRow,self.edgeSegDict_outerEast(1));
                   self.dl_outer(self.xHighRow,self.edgeSegDict_outerEast(end)), ...
                   self.dl_outer(self.yHighRow,self.edgeSegDict_outerEast(end))];

        end

        function val = get.edgeCoords_outerNorth(self)

            val = [self.dl_outer(self.xLowRow,self.edgeSegDict_outerNorth(1)), ...
                   self.dl_outer(self.yLowRow,self.edgeSegDict_outerNorth(1));
                   self.dl_outer(self.xHighRow,self.edgeSegDict_outerNorth(end)), ...
                   self.dl_outer(self.yHighRow,self.edgeSegDict_outerNorth(end))];

        end

        function val = get.edgeCoords_outerWest(self)

            val = [self.dl_outer(self.xLowRow,self.edgeSegDict_outerWest(1)), ...
                   self.dl_outer(self.yLowRow,self.edgeSegDict_outerWest(1));
                   self.dl_outer(self.xHighRow,self.edgeSegDict_outerWest(end)), ...
                   self.dl_outer(self.yHighRow,self.edgeSegDict_outerWest(end))];

        end

        function val = get.edgeCoords_outer(self)

            val(:,:,1) = self.edgeCoords_outerSouth;
            val(:,:,2) = self.edgeCoords_outerEast;
            val(:,:,3) = self.edgeCoords_outerNorth;
            val(:,:,4) = self.edgeCoords_outerWest;

        end

        function val = get.edgeCoords_inclusions(self)

            if self.nSeg_inclusions == 0, val = zeros(2,2,0);
            else,
                for i = 1:self.nSeg_inclusions
                    val(:,:,i) = [self.mat(self.xLowRow,self.edgeSegDict_inclusions{i}), ...
                                self.mat(self.yLowRow,self.edgeSegDict_inclusions{i});
                                self.mat(self.xHighRow,self.edgeSegDict_inclusions{i}), ...
                                self.mat(self.yHighRow,self.edgeSegDict_inclusions{i})];
                end
            end

        end

        function val = get.edgeCoords(self)

            val = cat(3,self.edgeCoords_outer,self.edgeCoords_inclusions);

        end


        % GETTERS: return segment data
        function val = segType(self,segID)

            val = self.mat(self.typeRow,segID);

        end

        function val = circleCenter(self,segID)
            
            val = [self.mat(self.centerxRow,segID), self.mat(self.centeryRow,segID)];

        end

        function val = circleRadius(self,segID)

            val = self.mat(self.radiusRow,segID);

        end

        function val = isLowerCircle(self,segID)

            val = mod(segID - self.nSeg_outer,4) == 1 | ...
                    mod(segID - self.nSeg_outer,4) == 2;
        
        end

        function val = isUpperCircle(self,segID)

            val = ~self.isLowerCircle(segID);
        
        end

        function val = squareEdgeOrientation(self,segID)
        % note: 1 == south, 2 == east, 3 == north, 4 == west
                
            val = mod(segID - self.nSeg_outer,4) + 1;
    
        end


        % ADD Y-LINES
        function self = add_yline(self,varargin)

			% collect values where yLines should be added
			if nargin == 1, yBar = mean(self.yLim);
			else, for i = 1:length(varargin), yBar(i) = varargin{i}; end, end
			
			% build new dl matrix
			self.mat = self.split_dl_outer(yBar);
			self.mat = self.increment_dl_outer(yBar);
			self.mat = self.increment_dl_inclusions(yBar);
			self.mat = self.set_dl_yLines;
	
        end
        
        function mat = split_dl_outer(self,yBar)

            % store dl_outer
            dl_outer = self.dl_outer;

			% split dl_outer columns
			for j = 1:length(yBar), 
				for i = [size(dl_outer,2):-1:1]
					if self.splitCol(yBar(j),dl_outer,i)
						
						% duplicate the i-th column
						dl_outer = [dl_outer(:,1:i), dl_outer(:,i:end)];
						
						% adjust y_start and y_stop
						dl_outer(5,i) = yBar(j);
						dl_outer(4,i+1) = yBar(j);
					end
				end
			end

            % store new dl matrix
            mat = [dl_outer self.dl_inclusions self.dl_yLines];

		end

		function mat = increment_dl_outer(self,yBar)

            % store dl_outer
            dl_outer = self.dl_outer;

            % increment dl_outer appropriately
			for j = 1:length(yBar)
				increment_cols = find(sum(dl_outer(4:5,:) > yBar(j),1) > 0);
				dl_outer(6,increment_cols) = dl_outer(6,increment_cols) + 1;
			end

            % store new dl matrix
            mat = [dl_outer self.dl_inclusions self.dl_yLines];

		end

		function mat = increment_dl_inclusions(self,yBar,dl_inclusions)

            % store dl_inclusions
            dl_inclusions = self.dl_inclusions;

			% increment region IDs
			for j = 1:length(yBar)
				increment_cols = find(sum(dl_inclusions(4:5,:) > yBar(j),1) > 0);
				dl_inclusions(7,increment_cols) = dl_inclusions(7,increment_cols) + 1;
			end

            % store new dl matrix
            mat = [self.dl_outer dl_inclusions self.dl_yLines];

		end

        function mat = set_dl_yLines(self)

			% find columns where y Line breaks up a vertical edge
			yChanges = find(self.dl_outer(4,1:end/2) ~= ...
										self.dl_outer(5,1:end/2));
			yChanges = yChanges(2:end);

            % build dl_ylines using these data
			dl_yLines =  ...
				[2 * ones(1,length(yChanges));				% line type
				self.xLim(1) * ones(1,length(yChanges));	% x_start
				self.xLim(2) * ones(1,length(yChanges));	% x_stop
				self.dl_outer(4,yChanges);					% y_start
				self.dl_outer(4,yChanges);					% y_stop
				self.dl_outer(6,yChanges);					% left region ID
				self.dl_outer(6,yChanges-1);				% right region ID
				zeros(3,length(yChanges))];			 	    % padding

            % store new dl matrix
            mat = [self.dl_outer, self.dl_inclusions, dl_yLines];

		end


        % INCLUSIONS ON/OFF
        function self = inclusionsON(self)

            % store dl_inclusions
            dl_inclusions = self.dl_inclusions;

            % number region IDs
            nInc = size(dl_inclusions,2) / 4;
            nInc = 1:nInc;
            nInc = repmat(nInc,4,1);
            nInc = reshape(nInc,1,[]);
            nInc = nInc + max(self.faceIDs_Omega_eps);

            % set region IDs to 1
            dl_inclusions(self.leftRow,:) = nInc;

            % store new dl matrix
            self.mat = [self.dl_outer dl_inclusions self.dl_yLines];

        end

        function self = inclusionsOFF(self)

            % store dl_inclusions
            dl_inclusions = self.dl_inclusions;

            % set left region IDs to 0
            dl_inclusions(self.leftRow,:) = 0;

            % store new dl matrix
            self.mat = [self.dl_outer dl_inclusions self.dl_yLines];

        end

    end

    methods (Static)

		function val = splitCol(yBar,dl_outer,i)

			val = (yBar > dl_outer(4,i) & yBar < dl_outer(5,i)) | ...
				(yBar > dl_outer(5,i) & yBar < dl_outer(4,i));
		
		end

    end

end