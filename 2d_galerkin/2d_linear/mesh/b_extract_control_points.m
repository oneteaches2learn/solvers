% Define the input and output file names
inputFile = 'trimmed.svg';
outputFile = 'controlpoints.txt';

% Read the content of the SVG file
fileContent = fileread(inputFile);

% Replace every 'M' with nothing
modifiedContent = regexprep(fileContent, 'M', '');

% Replace every 'L' with a newline followed by 'L'
modifiedContent = regexprep(fileContent, 'L', sprintf('\n'));

% Write the modified content to the new SVG file
fid = fopen(outputFile, 'w');
if fid == -1
    error('Cannot open file for writing: %s', outputFile);
end
fwrite(fid, modifiedContent);
fclose(fid);

disp(['Modified file saved as ', outputFile]);