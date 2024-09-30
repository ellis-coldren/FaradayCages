%create .poly test
filename = 'C:\Users\cathe\Documents\MATLAB\tetgen1.6.0\cubetest.poly';
fileID = fopen(filename, 'w');
% Cube -- boundary
bmax = 3;
bmin = -3;
boundary_vert = [bmin bmin bmin;
                bmin bmin bmax;
                bmin bmax bmax;
                bmin bmax bmin;
                bmax bmax bmax;
                bmax bmax bmin;
                bmax bmin bmin;
                bmax bmin bmax];


n = 1;  % resolution (number of subdivisions)
r = 1.5;
% -------------------- Circle ----------------------- %
[x_shape, y_shape, z_shape] = sphere(n);
c = [r*x_shape(:), r*y_shape(:), r*z_shape(:)];
c = unique(c,'rows');


% --------------------------------------------------- %

cage_points = [];
cage_vert_indices = [];
for i = 1:size(c,1)
    [wire_x, wire_y, wire_z] = sphere(5);
    wire_point = [1*wire_x(:)+c(i,1), 1*wire_y(:)+c(i,2), 1*wire_z(:)+c(i,3)];
    wire_point = unique(wire_point,'rows');
    wire_indices = convhull(wire_point);
    wire_indices = wire_indices + size(cage_points,1);
    column_of_i = i * ones(size(wire_indices, 1), 1);
    wire_indices = [wire_indices, column_of_i];
    cage_points = [cage_points; wire_point];
    cage_vert_indices = [cage_vert_indices;wire_indices];
end

points_per_wire = size(cage_points, 1) / size(c, 1);
figure;
trisurf(cage_vert_indices, cage_points(:,1), cage_points(:,2), cage_points(:,3), 'FaceColor', 'cyan', 'EdgeColor', 'black');
axis equal;


fprintf(fileID, '%d 3 0 1\n', size(boundary_vert, 1) + size(cage_points, 1));
for i = 1:size(boundary_vert, 1)
    fprintf(fileID, '%d %.6f %.6f %.6f %d\n', i, boundary_vert(i, 1), boundary_vert(i, 2), boundary_vert(i, 3), 1);
end
for i = 1:size(cage_points, 1)
    fprintf(fileID, '%d %.6f %.6f %.6f %d\n', i+size(boundary_vert, 1), cage_points(i, 1), cage_points(i, 2), cage_points(i, 3), 0);
end

boundary_faces = [
    1 2 3 4;  % face x = -1
    5 6 7 8;  % face x = 1
    1 2 8 7;  % face y = -1
    3 4 6 5;  % face y = 1
    1 4 6 7;  % face z = -1
    2 3 5 8;  % face z = 1
];

fprintf(fileID, '%d 1\n', size(boundary_faces, 1) + size(cage_vert_indices, 1));
disp(size(cage_vert_indices));

for i = 1:size(boundary_faces, 1)
    fprintf(fileID, '1 0 1 \n 4 %d %d %d %d #%d\n', boundary_faces(i, 1), boundary_faces(i, 2), boundary_faces(i, 3), boundary_faces(i, 4), i);
end
cage_vert_indices(:, 1:3) = cage_vert_indices(:, 1:3) + size(boundary_vert, 1);
count = size(boundary_faces, 1)+1;
for i = 1:size(cage_vert_indices, 1)
    fprintf(fileID, '1 0 %d\n 3 %d %d %d #%d\n', cage_vert_indices(i, 4)+1, cage_vert_indices(i, 1), cage_vert_indices(i, 2), cage_vert_indices(i, 3), i+size(boundary_faces,1));
end


fprintf(fileID, '%d\n', size(c,1));
for i = 1:size(c, 1)
    fprintf(fileID, '%d %d %d %d\n', i, c(i, 1), c(i, 2), c(i, 3));
end
fprintf(fileID, '1\n');
fprintf(fileID, '%d %d %d %d %d\n', 1, 8, 8, 8, 1);


fclose(fileID);