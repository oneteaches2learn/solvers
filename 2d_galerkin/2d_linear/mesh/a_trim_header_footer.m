% MATLAB Script to Modify SVG Content

% Define input and output file names
inputFile = 'raw.svg';
outputFile = 'trimmed.svg';

% Read the content of the input file
fileContent = fileread(inputFile);

% Remove all text up to and including 'd="'
startIdx = strfind(fileContent, 'd="');
if isempty(startIdx)
    error('The string ''d="'' was not found in the file.');
end
fileContent = fileContent(startIdx + length('d="'):end);

% Remove all text from 'sodipodi' to the end
endIdx = strfind(fileContent, 'sodipodi');
if ~isempty(endIdx)
    fileContent = fileContent(1:endIdx-1);
end

% Write the modified content to the output file
fid = fopen(outputFile, 'w');
if fid == -1
    error('Cannot open file for writing: %s', outputFile);
end
fwrite(fid, fileContent, 'char');
fclose(fid);

disp(['Modified file saved as ', outputFile]);