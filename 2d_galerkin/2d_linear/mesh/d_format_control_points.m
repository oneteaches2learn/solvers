% MATLAB Script to Convert 'modified_hand_controlpoints.svg' to 'points.geo'

% Define input and output file names
inputFile = 'filtered_controlpoints.txt';
outputFile = 'formatted_controlpoints.geo';

% Read the modified SVG file
fileContent = fileread(inputFile);

% Split the content into individual lines
lines = strsplit(fileContent, '\n');

% Open the output .geo file for writing
fid = fopen(outputFile, 'w');
if fid == -1
    error('Cannot open file for writing: %s', outputFile);
end

% Iterate over each line and format as Point(n) = {x, y, 0, 1.0};
for i = 1:length(lines)
    line = strtrim(lines{i});
    if isempty(line)
        continue; % Skip empty lines
    end
    coords = sscanf(line, '%f %f');
    if length(coords) >= 2
        fprintf(fid, 'Point(%d) = {%.4f, %.4f, 0, 1.0};\n', i, coords(1), coords(2));
    end
end

% Close the output file
fclose(fid);

disp(['Formatted points have been saved to ', outputFile]);