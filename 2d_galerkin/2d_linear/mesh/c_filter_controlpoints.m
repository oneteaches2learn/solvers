% Define the input file name
inputFile = 'controlpoints.txt';
outputFile = 'filtered_controlpoints.txt';

% Read the data into a matrix
data = readmatrix(inputFile);

% Check if data is not empty
if isempty(data)
    error('No valid data found in the input file.');
end

% Define the distance threshold (set by the user)
distanceThreshold = 5;  % Example threshold, adjust as needed

% Initialize a logical array to mark points for removal
toRemove = false(size(data, 1), 1);

% Iterate through the control points
i = 1;
while i <= size(data, 1)
    % Get the current control point
    currentPoint = data(i, :);
    
    % Check subsequent points
    for j = i+1:size(data, 1)
        nextPoint = data(j, :);
        % Calculate the Euclidean distance between the current and next point
        distance = sqrt((nextPoint(1) - currentPoint(1))^2 + (nextPoint(2) - currentPoint(2))^2);
        
        if distance <= distanceThreshold && j ~= size(data, 1)
            % Mark the next point for removal
            toRemove(j) = true;
        else
            % Move to the next unmarked control point
            break;
        end
    end
    
    % Move to the next unmarked control point
    i = find(~toRemove(i+1:end), 1) + i;
    if isempty(i)
        break;
    end
end

% Remove the marked points
filteredData = data(~toRemove, :);

% temporary: rescale the data
% 28534468085 is the factor that rescales the data from MP, PM, NS to the meter scale
filteredData = filteredData / 2853.4468085;

%{
% Normalize filtered data
normedData = norm(filteredData,2);
factor = max(normedData);
filteredData = filteredData/factor;
%}

% Display the filtered data
disp('Filtered Data:');
disp(filteredData);

% Optionally, save the filtered data to a new file
writematrix(filteredData, outputFile, 'Delimiter', 'tab');
disp(['Filtered data saved to ', outputFile]);