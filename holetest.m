% % % % n_points = 20;
% % % % 
% % % % % Radius of the sphere
% % % % r = 0.1;
% % % % 
% % % % % Generate spherical coordinates
% % % % theta = linspace(0, pi, n_points); % Polar angle
% % % % phi = linspace(0, 2*pi, n_points); % Azimuthal angle
% % % % 
% % % % % Create a grid for the angles
% % % % [Theta, Phi] = meshgrid(theta, phi);
% % % % 
% % % % % Convert spherical coordinates to Cartesian coordinates
% % % % x = r * sin(Theta) .* cos(Phi);
% % % % y = r * sin(Theta) .* sin(Phi);
% % % % z = r * cos(Theta);
% % % % 
% % % % xx = r * sin(theta) .* cos(phi);
% % % % yy = r * sin(theta) .* sin(phi);
% % % % zz = r * cos(theta);
% % % % 
% % % % vertices = [x(:), y(:), z(:)];
% % % % figure;
% % % % plot3(x(:), y(:), z(:), 'o');
% % % % xlabel('X');
% % % % ylabel('Y');
% % % % zlabel('Z');
% % % % axis equal;

n = 10;  % resolution (number of subdivisions)
[x, y, z] = sphere(n);
size(x)
size(y)
size(z)
% Reshape the coordinates into a list of points
points = [x(:), y(:), z(:)];

% Create the convex hull to get the triangulation of the sphere
size(points)
K = convhull(points);

% Plot the triangulated sphere
figure;
trisurf(K, points(:,1), points(:,2), points(:,3), 'FaceColor', 'cyan', 'EdgeColor', 'black');
axis equal;
