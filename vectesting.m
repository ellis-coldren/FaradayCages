n = 5;r = 0.01;
warning('off')
% -------------------- Reading .vec files ----------------------- %
vecFile = 'venndiagram.vec';

vecRead = xmlread(vecFile);

sketch_edges = vecRead.getElementsByTagName('edge');
x_shape = [];
y_shape = [];

for i = 0:sketch_edges.getLength-1
    i
    edge = sketch_edges.item(i);
    % the xywdense has format: (layer(int) x1, y1, w1, x2, y2, w2, ...);
    densesample = char(edge.getAttribute('curve'));
    sample_substring = regexp(densesample, '[+-]?\d*\.?\d+', 'match');
    sample_double = str2double(sample_substring);

    length_x = length(sample_double(2:3:end));
    length_y = length(sample_double(3:3:end));

    x_indices = round(linspace(1, length_x, n));
    y_indices = round(linspace(1, length_y, n));
    
    x_sample = sample_double(2:3:end);
    x_to_add = x_sample(x_indices);
    y_sample = sample_double(3:3:end);
    y_to_add = y_sample(y_indices);

    x_shape = [x_shape, x_sample];
    y_shape = [y_shape, y_sample];
end
min_x = min(x_shape);
max_x = max(x_shape);
min_y = min(y_shape);
max_y = max(y_shape);
x_scale = (3 - (-3))/(max_x - min_x);
y_scale = (3 - (-3))/(max_y - min_y);
x_shape = x_scale*x_shape;
y_shape = -y_scale*y_shape;


p = [x_shape', y_shape'];
% q = curvspace(p,n+1); % generates points that interpolate curve
% xx = q(:, 1); % x values of curve
% yy = q(:, 2); % y values of curve
xx = p(:, 1);
yy = p(:, 2);
c = xx + 1i*yy;
disp(c)

figure;
plot(real(c), imag(c), 'o', 'LineStyle', 'none');
xlabel('Real Part');
ylabel('Imaginary Part');
title('Plot of Complex Numbers');
grid on;


rr = r*ones(size(c)); % radii of cage vertices as vector
N = max(0, round(4+.5*log10(r)));
npts=3*N+2;
% cagePts_num = n*npts;
% xres = 5;
% yres = 5;



z = exp(pi*1i*(-50:50)'/50); % vector of circle coordinates for cage vertices
min_x = min(x_shape);
max_x = max(x_shape);
min_y = min(y_shape);
max_y = max(y_shape);
x = [min_x-2, max_x+2, max_x+2, min_x-2]; % boundary of region
y = [min_y-2, min_y-2, max_y+2, max_y+2]; % boundary of region

% https://www.mathworks.com/help/matlab/ref/polyshape.html %
% ---------------- setting up cage vertices --------------------- %

% connects vertices of x and y, adds hole where first cage vertex is (c(1) is the center of the cage point, + rr(1)*z draws a circle with radius rr(1) around it) %
pgon = polyshape({x, real(c(1)+rr(1)*z)}, {y, imag(c(1)+rr(1)*z)}); 
for j = 2:length(c)
    % adds another hole for each cage vertex %
    pgon = addboundary(pgon,real(c(j)+rr(j)*z),imag(c(j)+rr(j)*z));
end
figure;
plot(pgon);
 
% triangulates the region %
tr = triangulation(pgon);
figure;
triplot(tr);
points = tr.Points; %row number is vertex id
connectivity = tr.ConnectivityList; %row number is triangle id, rows contain vertex id

% creates finite element model so that we can create a mesh %
% https://www.mathworks.com/help/pde/ug/fegeometry.html

gm = fegeometry(tr);
gm = generateMesh(gm, 'Hmax', 0.5, 'Hmin', 0.001);
numelements = size(gm.Mesh.Elements, 2);
numnodes = size(gm.Mesh.Nodes, 2);
pdemesh(gm)
disp(numnodes);

neighborsMatrix = sparse(numnodes, numnodes);
