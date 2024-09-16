[status, cmdout] = system('tetgen -pqa1 C:\Users\cathe\Documents\MATLAB\tetgen1.6.0\cubetest.poly');

% Check if the command ran successfully
if status == 0
    disp('Command executed successfully:');
    disp(cmdout);
else
    disp('Error executing the command.');
    disp(cmdout);
end

% Specify the .ele file path
eleFile = 'C:\Users\cathe\Documents\MATLAB\tetgen1.6.0\cubetest.1.ele';

% Open the file
fid = fopen(eleFile, 'r');
if fid == -1
    error('File could not be opened. Check the file path.');
end

% Read the first line (header)
header = textscan(fid, '%d %d %d', 1);
numElements = header{1};  % Number of elements
nodesPerElement = header{2};  % Nodes per element (should be 4 for tetrahedra)

% Read the rest of the file (element IDs and their corresponding node indices)
formatSpec = '%d';  % Element ID
for i = 1:nodesPerElement
    formatSpec = strcat(formatSpec, ' %d');  % Node indices
end

% Read the elements
elements = textscan(fid, formatSpec, numElements);

% Close the file
fclose(fid);

% Extract element data (skip the first column if it's just element IDs)
elements = cell2mat(elements);
elementIDs = elements(:, 1);  % Element IDs
nodeIndices = elements(:, 2:end);  % Node indices for each element

nodeFile = 'C:\Users\cathe\Documents\MATLAB\tetgen1.6.0\cubetest.1.node';
fid = fopen(nodeFile, 'r');
if fid == -1
    error('File could not be opened. Check the file path.');
end
% Read the first line (header)
header = textscan(fid, '%d %d %d %d', 1);
numNodes = header{1};  % Number of elements
dim = header{2};
attrbNode = header{3};
boundarymarker = header{4};

% Read the rest of the file (element IDs and their corresponding node indices)
formatSpec = '%f';  % Element ID
for i = 1:(dim + attrbNode + boundarymarker)
    formatSpec = strcat(formatSpec, ' %f');  % Node indices
end

% Read the elements
nodes = textscan(fid, formatSpec, numNodes);
nodes = cell2mat(nodes);
% Close the file
fclose(fid);

% Visualize the tetrahedral mesh (assuming 3D)
tetramesh(nodeIndices, nodes(:, 2:4), 'FaceAlpha', 0.1, 'FaceColor', 'blue');  % nodeIndices are from the .ele file
title('Tetgen 3D Mesh');
xlabel('X');
ylabel('Y');
zlabel('Z');
axis equal;