% MATLAB Script to Replace '1.0};' with 'lc};' in 'mymodel2.geo'

% Define the input and output file names
inputFile = 'formatted_controlpoints.geo';
outputFile = 'formatted_controlpoints.geo';

% Read the content of the input file
fileContent = fileread(inputFile);

% Split the content into individual lines
lines = strsplit(fileContent, '\n');

% Initialize a cell array to hold modified lines
modifiedLines = cell(size(lines));

% Iterate over each line to perform the replacement
for i = 1:length(lines)
    line = lines{i};
    trimmedLine = strtrim(line);
    if startsWith(trimmedLine, 'Point(')
        % Replace '1.0};' with 'lc};'
        modifiedLine = regexprep(line, '1\.0};$', 'lc};');
        modifiedLines{i} = modifiedLine;
    else
        % Keep the line unchanged
        modifiedLines{i} = line;
    end
end

% Combine the modified lines back into a single string
modifiedContent = strjoin(modifiedLines, '\n');

% Write the modified content to the output file
fid = fopen(outputFile, 'w');
if fid == -1
    error('Cannot open file for writing: %s', outputFile);
end
fwrite(fid, modifiedContent, 'char');
fclose(fid);

disp(['Modified file saved as ', outputFile]);