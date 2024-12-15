classdef dlObj

    properties
        mat
        mat_reduced
    end

    properties (Dependent)
        nSeg_total
        nSeg_outer
        nSeg_inclusions
        nSeg_yLines
        nSeg_yLines_inc
        segIDs_outer
        segIDs_inclusions
        segIDs_yLines
        segIDs_yLines_inc
        dl_outer
        dl_inclusions
        dl_yLines
        dl_yLines_inc
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
        function self = dlObj(dl_mat)
            self.mat = self.set_dl(dl_mat);
            self.mat_reduced = self.set_dl_reduced;
        end


        % SETTERS
        function dl_mat = set_dl(self,dl_mat)

            % if all dl_mat cols correspond to lines, pad bottom of dl_mat
            nRows = size(dl_mat,1);
            if nRows ~= 10
                nCols = size(dl_mat,2);
                pad = zeros(10 - nRows,nCols);
                dl_mat = [dl_mat; pad];
            end
        end

        function self = set_dl_reduced(self)
        % note: this should be used only to initialize the dl_reduced matrix,
        % i.e. before the construction of any yLines. If this is used once
        % yLines exist, it will delete them and cause problems.
        
            % set reduced dl matrix
            self = self.reduce_dl_mat(self.dl_outer,self.dl_inclusions,[],[]);
    
        end


        % GETTERS: number of segments
        function val = get.nSeg_total(self)
        % primative getter, only should use data from dl_mat itself

            val = size(self.mat,2);

        end

        function val = get.nSeg_outer(self)
        % primative getter, only should use data from dl_mat itself
        
            val = sum(self.mat(self.rightRow,:) == 0);
         
        end

        function val = get.nSeg_yLines(self)
        % primative getter, only should use data from dl_mat itself
        %   note: yLines are defined by the following logic criteria:
        %       1. type == 2, corresponding to line segments
        %       2. yLow == yHigh, i.e. start and end at the same height
        %       3. yLow > yLowVal, i.e. above the bottom edge of the domain    
        %       4. yHigh < yHighVal, i.e. below the top edge of the domain
        %       5. left faceIDs are in the range of face IDs associated with outer edges
        %           (i.e. not within the range of those associated with inclusions)
            
            % get logic
            logic1 = self.mat(self.typeRow,:) == 2;
            logic2 = self.mat(self.yLowRow,:) == self.mat(self.yHighRow,:);
            logic3 = self.mat(self.yLowRow,:) > self.yLowVal;
            logic4 = self.mat(self.yHighRow,:) < self.yHighVal;
            logic5 = self.mat(self.leftRow,:) <= ...
                        max(max(self.mat(self.leftRow,1:self.nSeg_outer)), ...
                        max(self.mat(self.rightRow,1:self.nSeg_outer)));

            % combine logic
            logic_tot = logic1 & logic2 & logic3 & logic4 & logic5;

            % get total number of yLines
            val = sum(logic_tot);

        end

        function val = get.nSeg_yLines_inc(self)
        % primative getter, only should use data from dl_mat itself
        %   note: yLines are defined by the following logic criteria:
        %       1. type == 2, corresponding to line segments
        %       2. yLow == yHigh, i.e. start and end at the same height
        %       3. yLow > yLowVal, i.e. above the bottom edge of the domain    
        %       4. yHigh < yHighVal, i.e. below the top edge of the domain
        %       5. left faceIDs are outside the range of face IDs associated with outer edges
        %           (i.e. within the range of those associated with inclusions)
        %       6. right faceIDs are also outside the range of faceIDs associated with outer edges
        
            % get logic
            logic1 = self.mat(self.typeRow,:) == 2;
            logic2 = self.mat(self.yLowRow,:) == self.mat(self.yHighRow,:);
            logic3 = self.mat(self.yLowRow,:) > self.yLowVal;
            logic4 = self.mat(self.yHighRow,:) < self.yHighVal;
            logic5 = self.mat(self.leftRow,:) > ...
                        max(max(self.mat(self.leftRow,1:self.nSeg_outer)), ...
                        max(self.mat(self.rightRow,1:self.nSeg_outer)));
            logic6 = self.mat(self.rightRow,:) > ...
                        max(max(self.mat(self.leftRow,1:self.nSeg_outer)), ...
                        max(self.mat(self.rightRow,1:self.nSeg_outer)));

            % combine logic
            logic_tot = logic1 & logic2 & logic3 & logic4 & logic5 & logic6;

            % get total number of yLines
            val = sum(logic_tot);

        end

        function val = get.nSeg_inclusions(self)
        % primative getter, only should use data from dl_mat itself

            val = self.nSeg_total - self.nSeg_outer - self.nSeg_yLines - self.nSeg_yLines_inc;

        end

        % GETTERS: domain boundaries
        function val = get.xLowVal(self)
        % primative getter, only should use data from dl_mat itself
            val = min(self.mat(self.xLowRow:self.xHighRow,:),[],'all');
        end

        function val = get.xHighVal(self)
        % primative getter, only should use data from dl_mat itself
            val = max(self.mat(self.xLowRow:self.xHighRow,:),[],'all');
        end

        function val = get.yLowVal(self)
        % primative getter, only should use data from dl_mat itself
            val = min(self.mat(self.yLowRow:self.yHighRow,:),[],'all');
        end

        function val = get.yHighVal(self)
        % primative getter, only should use data from dl_mat itself
            val = max(self.mat(self.yLowRow:self.yHighRow,:),[],'all');
        end

        function val = get.xLim(self)
        % primative getter, only should use data from dl_mat itself
            val = [self.xLowVal self.xHighVal];
        end

        function val = get.yLim(self)
        % primative getter, only should use data from dl_mat itself
            val = [self.yLowVal self.yHighVal];
        end

        % GETTERS: segment IDs
        function val = get.segIDs_outer(self)
        % primative getter, only should use data from dl_mat itself
        
            val = find(self.mat(self.rightRow,:) == 0);
         
        end

        function val = get.segIDs_yLines(self)
        % primative getter, only should use data from dl_mat itself
        %   note: yLines are defined by the following logic criteria:
        %       1. type == 2, corresponding to line segments
        %       2. yLow == yHigh, i.e. start and end at the same height
        %       3. yLow > yLowVal, i.e. above the bottom edge of the domain    
        %       4. yHigh < yHighVal, i.e. below the top edge of the domain
        %       5. left faceIDs are in the range of face IDs associated with outer edges
        %           (i.e. not within the range of those associated with inclusions)

            % get logic
            logic1 = self.mat(self.typeRow,:) == 2;
            logic2 = self.mat(self.yLowRow,:) == self.mat(self.yHighRow,:);
            logic3 = self.mat(self.yLowRow,:) > self.yLowVal;
            logic4 = self.mat(self.yHighRow,:) < self.yHighVal;
            logic5 = self.mat(self.leftRow,:) <= ...
                        max(max(self.mat(self.leftRow,1:self.nSeg_outer)), ...
                        max(self.mat(self.rightRow,1:self.nSeg_outer)));

            % combine logic
            logic_tot = logic1 & logic2 & logic3 & logic4 & logic5;

            % get total number of yLines
            val = find(logic_tot);

        end

        function val = get.segIDs_yLines_inc(self)
        % primative getter, only should use data from dl_mat itself
        %   note: yLines_inc are defined by the following logic criteria:
        %       1. type == 2, corresponding to line segments
        %       2. yLow == yHigh, i.e. start and end at the same height
        %       3. yLow > yLowVal, i.e. above the bottom edge of the domain    
        %       4. yHigh < yHighVal, i.e. below the top edge of the domain
        %       5. left faceIDs are outside the range of face IDs associated with outer edges
        %           (i.e. within the range of those associated with inclusions)
        %       6. right faceIDs are also outside the range of faceIDs associated with outer edges

            % get logic
            logic1 = self.mat(self.typeRow,:) == 2;
            logic2 = self.mat(self.yLowRow,:) == self.mat(self.yHighRow,:);
            logic3 = self.mat(self.yLowRow,:) > self.yLowVal;
            logic4 = self.mat(self.yHighRow,:) < self.yHighVal;
            logic5 = self.mat(self.leftRow,:) > ...
                        max(max(self.mat(self.leftRow,1:self.nSeg_outer)), ...
                        max(self.mat(self.rightRow,1:self.nSeg_outer)));
            logic6 = self.mat(self.rightRow,:) > ...
                        max(max(self.mat(self.leftRow,1:self.nSeg_outer)), ...
                        max(self.mat(self.rightRow,1:self.nSeg_outer)));

            % combine logic
            logic_tot = logic1 & logic2 & logic3 & logic4 & logic5 & logic6;

            % get total number of yLines
            val = find(logic_tot);

        end

        function val = get.segIDs_inclusions(self)
        % depends on segIDs_outer, segIDs_yLines, segIDs_yLines_inc, each of
        % which are primative
        
            segIDs_total = 1:self.nSeg_total;
            temp = unique([self.segIDs_outer,self.segIDs_yLines,self.segIDs_yLines_inc]);
            val = setdiff(segIDs_total,temp);

        end


        % GETTERS: matrix partitions
        function val = get.dl_outer(self)
        % depends on segIDs_outer, which is primative

            val = self.mat(:,self.segIDs_outer);

        end

        function val = get.dl_inclusions(self)
        % depends on segIDs_inclusions

            val = self.mat(:,self.segIDs_inclusions);

        end

        function val = get.dl_yLines(self)
        % depends on segIDs_yLines, which is primative

            val = self.mat(:,self.segIDs_yLines);

        end

        function val = get.dl_yLines_inc(self)
        % depends on segIDs_yLines_inc, which is primative

            val = self.mat(:,self.segIDs_yLines_inc);

        end


        % GETTERS: segment-edge dictionary
        function val = get.segEdgeDict(self)
        % map between segment IDs and edge IDs
        
            val = [];
            for i = 1:length(self.edgeSegDict)
                
                val = [val i * ones(1,length(self.edgeSegDict{i}))];
                    
            end

        end


        % GETTERS: edge-segment dictionaries
        function val = get.edgeSegDict_outer(self)
        % map between edge IDs and segment IDs for outer edges/segments

            val{1} = self.edgeSegDict_outerSouth;
            val{2} = self.edgeSegDict_outerEast;
            val{3} = self.edgeSegDict_outerNorth;
            val{4} = self.edgeSegDict_outerWest;

        end

        function val = get.edgeSegDict_inclusions(self)
        % map between edge IDs and segment IDs for inclusion edges/segments

            % store data
            dl_inc = self.dl_inclusions;
            indices = 1:size(dl_inc,2);

            if ~isempty(indices)
                % get logic
                swLogic = dl_inc(self.yLowRow,:) >= dl_inc(self.yHighRow,:) & ...
                            dl_inc(self.xLowRow,:) < dl_inc(self.xHighRow,:);
                seLogic = dl_inc(self.yLowRow,:) < dl_inc(self.yHighRow,:) & ...
                            dl_inc(self.xLowRow,:) <= dl_inc(self.xHighRow,:);
                neLogic = dl_inc(self.yLowRow,:) <= dl_inc(self.yHighRow,:) & ...
                            dl_inc(self.xLowRow,:) > dl_inc(self.xHighRow,:);
                nwLogic = dl_inc(self.yLowRow,:) > dl_inc(self.yHighRow,:) & ...
                            dl_inc(self.xLowRow,:) >= dl_inc(self.xHighRow,:);

                % get indices
                swIndices = indices .* swLogic;
                seIndices = indices .* seLogic;
                neIndices = indices .* neLogic;
                nwIndices = indices .* nwLogic;

                % combine as rows of matrix
                indMat = [swIndices; seIndices; neIndices; nwIndices];
                indMat = indMat + self.nSeg_outer * (indMat > 0);

                % convert into cell array
                % note: this code moves along a row of the matrix, storing entries
                % as a vector until a zero is encountered. When a zero is
                % encountered, the vector is dumped into the next cell of a cell
                % array, and the loop moves to the next row of the matrix.
                row = 1;
                edge = 1;
                entries = indMat(1,1);
                for col = 2:length(indices)
                    if indMat(row,col) == 0
                        val{edge} = entries;
                        edge = edge + 1;
                        row = mod(row,4)+1;
                        entries = indMat(row,col);
                    else
                        entries = [entries,indMat(row,col)];
                    end
                end
                val{edge} = entries;
            
            else
                val = [];
            end

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
                edgeSegDict_inclusions = self.edgeSegDict_inclusions;
                for i = 1:length(edgeSegDict_inclusions)
                    entry = edgeSegDict_inclusions{i};
                    val(:,:,i) = [self.mat(self.xLowRow,entry(1)), ...
                                self.mat(self.yLowRow,entry(1));
                                self.mat(self.xHighRow,entry(end)), ...
                                self.mat(self.yHighRow,entry(end))];
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

        function val = isLowerCircle(self,edgeID)

            % for now, presume there are 4 outer edges
            val = mod(edgeID - 4,4) == 1 | ...
                    mod(edgeID - 4,4) == 2;
        
        end

        function val = isUpperCircle(self,segID)

            val = ~self.isLowerCircle(segID);
        
        end

        function val = squareEdgeOrientation(self,segID)
        % note: 1 == south, 2 == east, 3 == north, 4 == west
            
            % old code
            %val = mod(segID - self.nSeg_outer,4) + 1;

            % new code
            val = mod(segID - 4 - 1,4) + 1;
    
        end


        % ADD YLINE
        function self = add_yline(self,varargin)

			% collect values where yLines should be added
			if nargin == 1, yBar = mean(self.yLim);
			else, for i = 1:length(varargin), yBar(i) = varargin{i}; end, end
			
            % collect components
            dl_outer = self.dl_outer;
            dl_inc = self.dl_inclusions;
            dl_yLine = self.dl_yLines;
            yLine_height = yBar;

            % split outer boundary
            dl_outer = self.split_intersecting_columns(dl_outer,yLine_height);
            dl_outer = self.update_outer_faces(dl_outer,yLine_height);

            % split inclusions 
            dl_inc = self.split_intersecting_columns(dl_inc,yLine_height);
            dl_inc = self.update_inclusion_left_faces(dl_inc,yLine_height);
            dl_inc = self.update_inclusion_right_faces(dl_inc,yLine_height);

            % add yline
            [dl_yLine,dl_yLineInc] = self.add_dl_yLine(dl_outer,dl_inc,yLine_height);

            % reassemble dl matrix
            dl_mat = [dl_outer,dl_inc,dl_yLine,dl_yLineInc];

            % create reduced dl matrix
            dl_mat2 = self.reduce_dl_mat(dl_outer,dl_inc,dl_yLine,dl_yLineInc);

            % store results
            self.mat = dl_mat;
            self.mat_reduced = dl_mat2;

        end

        function intData = get_yLine_intersect_coords(self,dl_mat,yLine_height)
        % get_yLine_intersect_coords returns x and y coordinates of intersections

            % get logic 
            intCols_circ = dlObj.getCircularCols(dl_mat);
            intCols_line = dlObj.getLinearCols(dl_mat);

            % get intersection coordinates for circle
            intData.xCirc = self.get_circle_yLine_intersections(dl_mat,yLine_height);
            intData.yCirc = yLine_height * intCols_circ;

            % get intersection coordinates for square
            intData.xSqr = self.get_square_yLine_intersections(dl_mat,yLine_height);
            intData.ySqr = yLine_height * intCols_line;

            % combine into all shapes
            intData.xInt = intData.xCirc + intData.xSqr;
            intData.yInt = intData.yCirc + intData.ySqr;

        end

        function touchData = get_yLine_touch_coords(self,dl_mat,yLine_height)

            % get logic 
            touchingStart = dlObj.getTouchingStartCols(dl_mat,yLine_height);
            touchingEnd = dlObj.getTouchingEndCols(dl_mat,yLine_height);

            % replace zeros in touchData with NaN
            touchingStart = double(touchingStart);
            touchingEnd = double(touchingEnd);
            touchingStart(touchingStart == 0) = NaN;
            touchingEnd(touchingEnd == 0) = NaN;

            % get touch coordinates for circle
            touchData.xStarts = unique(dl_mat(dlObj.xLowRow,:) .* touchingStart);
            touchData.xEnds = unique(dl_mat(dlObj.xHighRow,:) .* touchingEnd);
            touchData.x = unique([touchData.xStarts, touchData.xEnds]);

            % omit nan
            touchData.xStarts = touchData.xStarts(~isnan(touchData.xStarts));
            touchData.xEnds = touchData.xEnds(~isnan(touchData.xEnds));
            touchData.x = touchData.x(~isnan(touchData.x));

        end

        function xCirc = get_circle_yLine_intersections(self,dl_mat,yLine_height)
        % get_circle_yLine_intersections returns x-coordinates of intersections

            % get logic 
            [intCols,intCols_asc,intCols_desc] = dlObj.getAllIntersectingCols(dl_mat,yLine_height);
            circCols = dlObj.getCircularCols(dl_mat);

            % get data for computing intersections
            x_c = dl_mat(dlObj.centerxRow,:) .* (intCols & circCols);
            y_c = dl_mat(dlObj.centeryRow,:) .* (intCols & circCols);
            r   = dl_mat(dlObj.radiusRow,:)  .* (intCols & circCols);
            h   = yLine_height * (intCols & circCols);

            % compute and store intersections
            xCirc_asc = (x_c + sqrt(r.^2 - (h - y_c).^2)) .* intCols_asc;
            xCirc_desc = (x_c - sqrt(r.^2 - (h - y_c).^2)) .* intCols_desc;
            xCirc = xCirc_asc + xCirc_desc;

        end

        function xSqr = get_square_yLine_intersections(self,dl_mat,yLine_height)
        % get_circle_yLine_intersections returns x-coordinates of intersections

            % get logic 
            intCols = dlObj.getAllIntersectingCols(dl_mat,yLine_height);
            lineCols = dlObj.getLinearCols(dl_mat);

            % square edge intersections
            xSqr = dl_mat(dlObj.xHighRow,:) .* (intCols & lineCols);

        end

        function dl_mat = split_intersecting_columns(self,dl_mat,yLine_height)
        % split_intersecting_columns duplicates columns intersecting yLine and updates the coordinates

            % get logic 
            intCols = dlObj.getAllIntersectingCols(dl_mat,yLine_height);

            % get coordinates 
            intData = self.get_yLine_intersect_coords(dl_mat,yLine_height);
            
            % duplicate intersecting columns
            duped_indices = sort([[1:size(dl_mat,2)], find(intCols)]);
            dl_mat = dl_mat(:,duped_indices);

            % duplicate intersection data appropriately
            xInt_duped = intData.xInt(:,duped_indices);
            yInt_duped = intData.yInt(:,duped_indices);

            % organize into start/end points
            newEnds = [diff(duped_indices) == 0,0];
            newStarts = [0,newEnds(:,1:end-1)];
            x_newEnd   = xInt_duped .* newEnds;
            x_newStart = xInt_duped .* newStarts;
            y_newEnd   = yInt_duped .* newEnds;
            y_newStart = yInt_duped .* newStarts;
            x = [x_newStart; x_newEnd];
            y = [y_newStart; y_newEnd];

            % zero out entries that should be overwritten, and then add in the new data
            dl_mat(dlObj.xLowRow:dlObj.xHighRow,:) = ...
                    dl_mat(dlObj.xLowRow:dlObj.xHighRow,:) .* ~x + x;
            dl_mat(dlObj.yLowRow:dlObj.yHighRow,:) = ...
                    dl_mat(dlObj.yLowRow:dlObj.yHighRow,:) .* ~y + y;

        end

        function dl_mat = update_outer_faces(self,dl_mat,yLine_height)

            % get logic 
            strictlyAbove = dlObj.getStrictlyAboveCols(dl_mat,yLine_height);
            touchingAbove = dlObj.getTouchingAboveCols(dl_mat,yLine_height);

            % default is to not increment
            faceIncrem = strictlyAbove | touchingAbove;

            % apply increments
            dl_mat(dlObj.leftRow,:) = dl_mat(dlObj.leftRow,:) + faceIncrem;

        end

        function dl_mat = update_inclusion_left_faces(self,dl_mat,yLine_height)

            % get logic 
            incTouch = dlObj.getInclusionsTouchingYLine(dl_mat,yLine_height);
            strictlyBelow = dlObj.getStrictlyBelowCols(dl_mat,yLine_height);
            touchingBelow = dlObj.getTouchingBelowCols(dl_mat,yLine_height);
            touchingAbove = dlObj.getTouchingAboveCols(dl_mat,yLine_height);
            strictlyAbove = dlObj.getStrictlyAboveCols(dl_mat,yLine_height);
            incEnum = dlObj.enumerateInclusions(dl_mat);
            incAbove = dlObj.getStrictlyAboveYLineCols(dl_mat,yLine_height);

            % columns corresponding to inclusions that straddle the y-line
            touchBelow = (incTouch & (touchingBelow | strictlyBelow));
            touchAbove = (incTouch & (touchingAbove | strictlyAbove));

            % number of first inclusion to touch y-line
            incNum = incEnum(...
                min([find(touchBelow,1,'first'),find(touchAbove,1,'first')]));
            
            % increment all inclusions by 1
            faceIncrem = ones(1,size(dl_mat,2));
            
            % increment those inclusions that are touching the y-line
            if sum(touchBelow) > 0 && sum(touchAbove) > 0
                faceIncrem = faceIncrem + (incEnum - incNum) .* touchBelow;
                faceIncrem = faceIncrem + (incEnum - incNum + 1) .* touchAbove;
            end

            % increment those inclusions that are strictly above the y-line
            faceIncrem = faceIncrem + (max(faceIncrem) - 1) .* incAbove;

            % apply increments
            dl_mat(dlObj.leftRow,:) = dl_mat(dlObj.leftRow,:) + faceIncrem;

        end

        function dl_mat = update_inclusion_right_faces(self,dl_mat,yLine_height)

            % get logic 
            strictlyAbove = dlObj.getStrictlyAboveCols(dl_mat,yLine_height);
            touchingAbove = dlObj.getTouchingAboveCols(dl_mat,yLine_height);

            % default is to not increment
            faceIncrem = strictlyAbove | touchingAbove;

            % apply increments
            dl_mat(dlObj.rightRow,:) = dl_mat(dlObj.rightRow,:) + faceIncrem;

        end

        function dl_mat = reduce_dl_mat(self,dl_outer,dl_inc,dl_yLine,dl_yLineInc)

            % set faceID inside inclusionts to 0
            dl_inc(dlObj.leftRow,:) = 0;

            % remove columns corresponding to line segments inside inclusions
            dl_mat = [dl_outer,dl_inc,dl_yLine];

        end

        function [dl_yLine,dl_yLineInc] = add_dl_yLine(self,dl_outer,dl_inc,yLine_height)
        % split_intersecting_columns duplicates columns intersecting yLine and updates the coordinates

            % get logic
            dl_input = [dl_outer,dl_inc];
            touchingAbove = dlObj.getTouchingAboveCols(dl_input,yLine_height);
            touchingBelow = dlObj.getTouchingBelowCols(dl_input,yLine_height);

            % get points where yLine touches shapes
            touchData = self.get_yLine_touch_coords(dl_input,yLine_height);

            % build yLine
            colNum = ones(1,length(touchData.x)-1);
            dl_mat =  ...
            [2 * colNum;                         % line type
                touchData.x(1:end-1);               % x_start
                touchData.x(2:end);                 % x_stop
                yLine_height * colNum;			    % y_start
                yLine_height * colNum;			    % y_stop
                0 * colNum;			                % left region ID
                0 * colNum;				            % right region ID
                0 * repmat(colNum,3,1)];	        % padding

            % if yLine is tangent to circles
            if sum(touchingAbove(size(dl_outer,2)+1:end)) == 0 || ...
                    sum(touchingBelow(size(dl_outer,2)+1:end)) == 0
                dl_input = dl_outer;

            % else, yLine is not tangent to circles
            else
                dl_input = [dl_outer,dl_inc];
            end

            % recalculate logic
            touchingBelow = dlObj.getTouchingBelowCols(dl_input,yLine_height);
            touchingAbove = dlObj.getTouchingAboveCols(dl_input,yLine_height);

            % get unique face IDs
            lowerFaceIDs = unique(touchingBelow .* dl_input(dlObj.leftRow,:));
            lowerFaceIDs = lowerFaceIDs(lowerFaceIDs ~= 0);
            upperFaceIDs = unique(touchingAbove .* dl_input(dlObj.leftRow,:));
            upperFaceIDs = upperFaceIDs(upperFaceIDs ~= 0);

            % reshape into appropriate array
            if length(lowerFaceIDs) > 1
                temp = [lowerFaceIDs(1) * ones(1,length(lowerFaceIDs)-1);
                        lowerFaceIDs([2:end])];
                lowerFaceIDs = reshape(temp,1,[]);
                lowerFaceIDs = [lowerFaceIDs,lowerFaceIDs(1)];
            else
                lowerFaceIDs = lowerFaceIDs(1) * ones(1,size(dl_mat,2));
            end

            if length(upperFaceIDs) > 1
                temp = [upperFaceIDs(1) * ones(1,length(upperFaceIDs)-1);
                        upperFaceIDs([2:end])];
                upperFaceIDs = reshape(temp,1,[]);
                upperFaceIDs = [upperFaceIDs,upperFaceIDs(1)];
            else
                upperFaceIDs = upperFaceIDs(1) * ones(1,size(dl_mat,2));
            end

            % assign face IDs
            dl_mat(dlObj.rightRow,:) = lowerFaceIDs;
            dl_mat(dlObj.leftRow,:) = upperFaceIDs;

            % move cols corresponding to line segs inside inclusions
            maxInteriorFaceID = max( ...
                max(dl_outer(dlObj.leftRow,:)),max(dl_outer(dlObj.rightRow,:)));
            moveCols = find(dl_mat(dlObj.leftRow,:) > maxInteriorFaceID);
            stayCols = setdiff(1:size(dl_mat,2),moveCols);

            % check for duplicates
            duplicates1 = intersect(dl_mat(1:5,:)',dl_inc(1:5,:)','rows');
            duplicates1 = duplicates1';

            % reverse x-coords and check for duplicates
            duplicates2 = dl_mat(1:5,:);
            duplicates2([2,3],:) = duplicates2([3,2],:);
            duplicates2 = intersect(duplicates2',dl_inc(1:5,:)','rows');
            duplicates2 = duplicates2';
            duplicates2([2,3],:) = duplicates2([3,2],:);

            % combine results
            duplicates = [duplicates1,duplicates2];

            % find indices of columns to delete
            deleteCols = zeros(1,size(duplicates,2));
            for i = 1:size(duplicates,2)
                deleteCols(i) = find(prod(dl_mat(1:5,:) == duplicates(:,i),1));
            end

            % save results
            dl_yLine = dl_mat(:,stayCols);
            dl_yLineInc = dl_mat(:,moveCols);

            % delete the columns
            dl_yLine(:,deleteCols) = [];

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

        function intCols_asc = getAscendingIntersectingCols(dl_mat,yLine_height)
            % get cols corresponding to ascending segments that intersect the y-line
            logic1 = dl_mat(dlObj.yLowRow, :) < yLine_height;
            logic2 = dl_mat(dlObj.yHighRow, :) > yLine_height;
            intCols_asc = logic1 & logic2;
        end

        function intCols_desc = getDescendingIntersectingCols(dl_mat,yLine_height)
            % get cols corresponding to descending segments that intersect the y-line
            logic1 = dl_mat(dlObj.yLowRow, :) > yLine_height;
            logic2 = dl_mat(dlObj.yHighRow, :) < yLine_height;
            intCols_desc = logic1 & logic2;
        end

        function [intCols,intCols_asc,intCols_desc] = getAllIntersectingCols(dl_mat,yLine_height)
            % get all cols intersecting yLine
            intCols_asc = dlObj.getAscendingIntersectingCols(dl_mat,yLine_height);
            intCols_desc = dlObj.getDescendingIntersectingCols(dl_mat,yLine_height);
            intCols = intCols_asc | intCols_desc;
        end

        function circCols = getCircularCols(dl_mat)
            % get cols corresponding to circular segments
            circCols = dl_mat(dlObj.typeRow, :) == 1;
        end

        function lineCols = getLinearCols(dl_mat)
            % get cols corresponding to linear segments
            lineCols = dl_mat(dlObj.typeRow, :) == 2;
        end

        function intCols_circ = getIntersectingCircularCols(dl_mat,yLine_height)
            % get cols intersecting yLine by shape (circular)
            intCols = dlObj.getAllIntersectingCols(dl_mat,yLine_height);
            circCols = dlObj.getCircularCols(dl_mat);
            intCols_circ = intCols & circCols;
        end

        function intCols_line = getIntersectingLinearCols(dl_mat,yLine_height)
            % get cols intersecting yLine by shape (linear)
            intCols = dlObj.getAllIntersectingCols(dl_mat,yLine_height);
            lineCols = dlObj.getLinearCols(dl_mat);
            intCols_line = intCols & lineCols;
        end

        function strictlyAbove = getStrictlyAboveCols(dl_mat,yLine_height)
            % get cols strictly above yLine
            logic1 = dl_mat(dlObj.yLowRow, :) > yLine_height;
            logic2 = dl_mat(dlObj.yHighRow, :) > yLine_height;
            strictlyAbove = logic1 & logic2;
        end
        
        function touchingAbove = getTouchingAboveCols(dl_mat,yLine_height)
            % get cols above but touching yLine
            logic1 = dl_mat(dlObj.yLowRow,:) == yLine_height;
            logic2 = dl_mat(dlObj.yHighRow,:) > yLine_height;
            logic3 = dl_mat(dlObj.yHighRow,:) == yLine_height;
            logic4 = dl_mat(dlObj.yLowRow,:) > yLine_height;
            touchingAbove = (logic1 & logic2) | (logic3 & logic4);
        end
        
        function touchingBelow = getTouchingBelowCols(dl_mat,yLine_height)
            % get cols below but touching yLine
            logic1 = dl_mat(dlObj.yLowRow,:) == yLine_height;
            logic2 = dl_mat(dlObj.yHighRow,:) < yLine_height;
            logic3 = dl_mat(dlObj.yHighRow,:) == yLine_height;
            logic4 = dl_mat(dlObj.yLowRow,:) < yLine_height;
            touchingBelow = (logic1 & logic2) | (logic3 & logic4);
        end
        
        function strictlyBelow = getStrictlyBelowCols(dl_mat,yLine_height)
            % get cols strictly below yLine
            logic1 = dl_mat(dlObj.yLowRow,:) < yLine_height;
            logic2 = dl_mat(dlObj.yHighRow,:) < yLine_height;
            strictlyBelow = logic1 & logic2;
        end
        
        function touchingStart = getTouchingStartCols(dl_mat,yLine_height)
            % get cols with start touching yLine
            touchingStart = dl_mat(dlObj.yLowRow,:) == yLine_height;
        end
        
        function touchingEnd = getTouchingEndCols(dl_mat,yLine_height)
            % get cols with end touching yLine
            touchingEnd = dl_mat(dlObj.yHighRow,:) == yLine_height;
        end
        
        function newInc = detectFirstColOfEachInclusion(dl_mat)
            % detect first col of each inclusion
            newInc = [1, dl_mat(dlObj.xHighRow, 1:end-1) ~= dl_mat(dlObj.xLowRow, 2:end)];
        end
        
        function incEnum = enumerateInclusions(dl_mat)
            % enumerate inclusions
            newInc = dlObj.detectFirstColOfEachInclusion(dl_mat);
            incEnum = zeros(1, length(newInc));
            incEnum(1) = 1;
            for i = 2:length(newInc), incEnum(i) = incEnum(i-1) + newInc(i); end
        end
        
        function incTouch = getInclusionsTouchingYLine(dl_mat,yLine_height)
            % get numbers for inclusions touching the y-line
            incEnum = dlObj.enumerateInclusions(dl_mat);
            touchingBelow = dlObj.getTouchingBelowCols(dl_mat,yLine_height);
            touchingAbove = dlObj.getTouchingAboveCols(dl_mat,yLine_height);
            incNums = unique(incEnum(find(incEnum .* (touchingBelow | touchingAbove))));
            incTouch = zeros(1, length(incEnum));
            for i = 1:length(incNums)
                incTouch = incTouch + (incEnum == incNums(i));
            end
        end
        
        function incBelow = getStrictlyBelowYLineCols(dl_mat,yLine_height)
            % get columns corr. to inc's strictly below y-line
            incTouch = dlObj.getInclusionsTouchingYLine(dl_mat,yLine_height);
            strictlyBelow = dlObj.getStrictlyBelowCols(dl_mat,yLine_height);
            incBelow = strictlyBelow .* ~incTouch;
        end
        
        function incAbove = getStrictlyAboveYLineCols(dl_mat,yLine_height)
            % get columns corr. to inc's strictly above y-line
            incTouch = dlObj.getInclusionsTouchingYLine(dl_mat,yLine_height);
            strictlyAbove = dlObj.getStrictlyAboveCols(dl_mat,yLine_height);
            incAbove = strictlyAbove .* ~incTouch;
        end

    end
end