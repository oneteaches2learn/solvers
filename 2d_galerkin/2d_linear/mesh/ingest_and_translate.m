% Define the input file name
inputFile = 'controlpoints.txt';

% Read the data into a matrix
data = readmatrix(inputFile);

% Check if data is not empty
if isempty(data)
    error('No valid data found in the input file.');
end

% Find the minimum value in each column
minX = min(data(:, 1));
minY = min(data(:, 2));

% Subtract the minimum value from each entry in the respective columns
data(:, 1) = data(:, 1) - minX;
data(:, 2) = data(:, 2) - minY;

% Display the normalized data
disp('Normalized Data:');
disp(data);

% Optionally, save the normalized data to a new file
outputFile = 'normalized_controlpoints.txt';
writematrix(data, outputFile, 'Delimiter', 'tab');
disp(['Normalized data saved to ', outputFile]);