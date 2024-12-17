% MATLAB Script to Convert 'points.geo' to 'mymodel2.geo' with Line Definitions

% Define input and output file names
inputFile = 'formatted_controlpoints.geo';
outputFile = 'domain.geo';

% Read the content of the points.geo file
fileContent = fileread(inputFile);

% Split the content into individual lines
lines = strsplit(fileContent, '\n');

% Extract Point lines
pointLines = lines(startsWith(strtrim(lines), 'Point('));

% Count the number of points
numPoints = length(pointLines);

% Open the output file for writing
fid = fopen(outputFile, 'w');
if fid == -1
    error('Cannot open file for writing: %s', outputFile);
end

% Write the original points to the new file
for i = 1:length(lines)
    fprintf(fid, '%s\n', lines{i});
end

% Add Line definitions
for i = 1:numPoints-1
    fprintf(fid, 'Line(%d) = {%d, %d};\n', i, i, i+1);
end
fprintf(fid, 'Line(%d) = {%d, %d};\n', numPoints, numPoints, 1);

% Add Curve Loop definition
fprintf(fid, '\nCurve Loop(%d) = {', numPoints + 1);
fprintf(fid, '%d, ', 1:numPoints);
fprintf(fid, '};\n');

% Add Plane Surface definition
fprintf(fid, 'Plane Surface(%d) = {%d};\n', numPoints + 2, numPoints + 1);

% Close the file
fclose(fid);

disp(['mymodel2.geo has been created with ', num2str(numPoints + 2), ' definitions.']);