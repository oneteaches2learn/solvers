% MATLAB Script to Delete Even Points and Connect Odd Points in 'mymodel2_modified.geo'

% Define input and output file names
inputFile = 'formatted_controlpoints.geo';
outputFile = 'domain.geo';

% Read the content of the input file
fileContent = fileread(inputFile);

% Split the content into individual lines
lines = strsplit(fileContent, '\n');

% Initialize cell arrays to hold different sections
headerLines = {};      % To store lines before Point definitions
pointLines = {};       % To store Point definitions
otherLines = {};       % To store lines after Point definitions

% Flags to determine the parsing state
inPointSection = false;

% Iterate through each line to separate header, points, and other content
for i = 1:length(lines)
    line = strtrim(lines{i});
    
    if startsWith(line, 'Point(')
        inPointSection = true;
        pointLines{end+1} = lines{i}; %#ok<AGROW>
    elseif inPointSection && isempty(line)
        inPointSection = false;
        otherLines{end+1} = lines{i}; %#ok<AGROW>
    elseif inPointSection
        pointLines{end+1} = lines{i}; %#ok<AGROW>
    else
        headerLines{end+1} = lines{i}; %#ok<AGROW>
    end
end

% Process Point definitions: Keep only odd-numbered points
keptPoints = {};
originalToNewMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
newPointIndex = 1;

for i = 1:length(pointLines)
    line = strtrim(pointLines{i});
    % Use regular expression to extract point number and coordinates
    tokens = regexp(line, 'Point\((\d+)\)\s*=\s*\{([^,]+),\s*([^,]+),\s*([^,]+),\s*([^}]+)\};', 'tokens');
    if ~isempty(tokens)
        pointNum = str2double(tokens{1}{1});
        x = tokens{1}{2};
        y = tokens{1}{3};
        z = tokens{1}{4};
        lc = tokens{1}{5};
        
        if mod(pointNum, 2) == 1  % Keep odd points
            % Map original point number to new point number
            originalToNewMap(pointNum) = newPointIndex;
            % Create new Point definition
            newPointLine = sprintf('Point(%d) = {%s, %s, %s, %s};', newPointIndex, x, y, z, lc);
            keptPoints{end+1} = newPointLine; %#ok<AGROW>
            newPointIndex = newPointIndex + 1;
        end
    end
end

% Total number of kept points
numKeptPoints = length(keptPoints);

% Generate Line definitions connecting the kept points sequentially
keptLines = {};
for i = 1:numKeptPoints-1
    keptLines{end+1} = sprintf('Line(%d) = {%d, %d};', i, i, i+1); %#ok<AGROW>
end
% Connect the last point back to the first to close the loop
keptLines{end+1} = sprintf('Line(%d) = {%d, %d};', numKeptPoints, numKeptPoints, 1); %#ok<AGROW>

% Generate Curve Loop definition
curveLoopId = numKeptPoints + 1;
curveLoopDef = sprintf('Curve Loop(%d) = {', curveLoopId);
for i = 1:numKeptPoints
    curveLoopDef = [curveLoopDef, sprintf('%d, ', i)]; %#ok<AGROW>
end
curveLoopDef(end-1:end) = '};';  % Replace the last comma and space with closing brace
keptLines{end+1} = curveLoopDef;

% Generate Plane Surface definition
planeSurfaceId = curveLoopId + 1;
planeSurfaceDef = sprintf('Plane Surface(%d) = {%d};', planeSurfaceId, curveLoopId);
keptLines{end+1} = planeSurfaceDef;

% Write the modified content to the output file
fid = fopen(outputFile, 'w');
if fid == -1
    error('Cannot open file for writing: %s', outputFile);
end

% Write header lines
for i = 1:length(headerLines)
    fprintf(fid, '%s\n', headerLines{i});
end

% Write kept Point definitions
for i = 1:length(keptPoints)
    fprintf(fid, '%s\n', keptPoints{i});
end

% Write Line definitions
for i = 1:length(keptLines)
    fprintf(fid, '%s\n', keptLines{i});
end

% Write other lines (if any)
for i = 1:length(otherLines)
    fprintf(fid, '%s\n', otherLines{i});
end

% Close the file
fclose(fid);

disp(['Modified file with odd points and connecting lines saved as ', outputFile]);