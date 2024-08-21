n = 15;r = 0.01;

t = linspace(0, 2*pi, 1000);

%Circle
% radius = 0.75;
% x_shape = radius*cos(t);
% y_shape = radius*sin(t);


%Heart
x_shape = 0.08*(16*sin(t).^3);
y_shape = 0.08*(13*cos(t) - 5*cos(2*t) - 2*cos(3*t) - cos(4*t));

%Diamond
% a = 0.7;
% x_shape = a*cos(t).^3;
% y_shape = a*sin(t).^3;

%Rose Curve
% a = 0.6;
% petals = 4;
% x_shape = a*cos(petals.*t).*cos(t);
% y_shape = a*cos(petals.*t).*sin(t);

%Hypotrochoid(sometimes playing around messes up, see https://www.desmos.com/calculator/r8bsgvsszo)
% fr = 1.6;
% rr = 0.2;
% d = 0.5;
% x_shape = 0.5*((fr-rr)*cos(t)+d*cos(((fr-rr)/rr)*t));
% y_shape = 0.5*((fr-rr)*sin(t) - d*sin(((fr-rr)/rr)*t));


p = [x_shape', y_shape'];

q = curvspace(p,n+1);
%disp(q);
xx = q(:, 1);
yy = q(:, 2);
c = xx + 1i*yy;
disp(c);
disp(size(c));

% c = exp(2i*pi*(1:n)/n);
rr = r*ones(size(c));
N = max(0, round(4+.5*log10(r)));
npts=3*N+2;
cagePts_num = n*npts;
xres = 5;
yres = 5;

% x = linspace(-3, 3, 10); y = linspace(-3, 3, 10);


z = exp(pi*1i*(-50:50)'/50);
% x = [0.25, 0.75, 0.75, 0.25];
% y = [0.6, 0.6, 1.2, 1.2];
x = [-3, 3, 3, -3];
y = [-3, -3, 3, 3];
disp(c);
pgon = polyshape({x, real(c(1)+rr(1)*z)}, {y, imag(c(1)+rr(1)*z)});
for j = 2:n
    pgon = addboundary(pgon,real(c(j)+rr(j)*z),imag(c(j)+rr(j)*z));
end
tr = triangulation(pgon);
points = tr.Points; %row number is vertex id
connectivity = tr.ConnectivityList; %row number is triangle id, rows contain vertex id
% neighbors = neighbors(tr, id); %row is triangle id, tows contain triangle ids
gm = fegeometry(tr);
gm = generateMesh(gm, Hmax = 0.2, Hmin= 0.001);
numelements = size(gm.Mesh.Elements, 2);
numnodes = size(gm.Mesh.Nodes, 2);
% midpoints of the edge are the last 3 rows of the gm.Mesh.Elements
% vertices of the element(triangle) are the first 3 rows of the
% gm.Mesh.Elements

% finding cage nodes:
cagePts = findNodes(gm.Mesh,"radius",[real(c(1)), imag(c(1))], rr(1));
for j = 1:n
    cagePts = [cagePts findNodes(gm.Mesh,"radius",[real(c(j)), imag(c(j))], rr(j))];
end

boundaryPts = findNodes(gm.Mesh, "region", "Edge",[1 2 3 4]);

neighborsMatrix = zeros(numnodes);
cotangents = zeros(numnodes);

for i = 1:numelements
  v1_id = gm.Mesh.Elements(1, i);
  v2_id = gm.Mesh.Elements(2, i);
  v3_id = gm.Mesh.Elements(3, i);

  v1 = gm.Mesh.Nodes(:, v1_id);
  v2 = gm.Mesh.Nodes(:, v2_id);
  v3 = gm.Mesh.Nodes(:, v3_id);
  
  neighborsMatrix(v1_id, v2_id) = 1;
  cotangents(v1_id, v2_id) = cotangents(v1_id, v2_id) + 0.5*cot(acos(dot(v1-v3, v2-v3)/(norm(v1-v3)*norm(v2-v3))));
  cotangents(v2_id, v1_id) = cotangents(v1_id, v2_id);
  neighborsMatrix(v2_id, v3_id) = 1;
  cotangents(v2_id, v3_id) = cotangents(v2_id, v3_id) + 0.5*cot(acos(dot(v2-v1, v3-v1)/(norm(v2-v1)*norm(v3-v1))));
  cotangents(v3_id, v2_id) = cotangents(v2_id, v3_id);
  neighborsMatrix(v3_id, v1_id) = 1;
  cotangents(v3_id, v1_id) = cotangents(v3_id, v1_id) + 0.5*cot(acos(dot(v3-v2, v1-v2)/(norm(v3-v2)*norm(v1-v2))));
  cotangents(v1_id, v3_id) = cotangents(v3_id, v1_id);
  neighborsMatrix(v2_id, v1_id) = 1;
  neighborsMatrix(v3_id, v2_id) = 1;
  neighborsMatrix(v1_id, v3_id) = 1;
end

neighborsMatrix = neighborsMatrix(any(neighborsMatrix, 2), :); %removing rows of the middle point nodes
neighborsMatrix = neighborsMatrix(:, any(neighborsMatrix, 1)); %removing columns of the middle point nodes

cotangents = cotangents(any(cotangents, 2), :); %removing rows of the middle point nodes
cotangents = cotangents(:, any(cotangents, 1)); %removing columns of the middle point nodes
size(neighborsMatrix)
numnodes = size(neighborsMatrix, 1);

colstodelete = cagePts > numnodes;
cagePts(colstodelete) = [];
colstodelete = boundaryPts > numnodes;
boundaryPts(colstodelete) = [];

% Laplacian_umbrella = zeros(numnodes);
Laplacian_cotangents = zeros(numnodes);
for i = 1:numnodes
    % normalize_value1 = sum(neighborsMatrix(i, :));
    normalize_value2 = sum(cotangents(i, :));
    % Laplacian_umbrella(i, :) = -neighborsMatrix(i, :)/normalize_value1;
    % Laplacian_umbrella(i:i) = 1;
    Laplacian_cotangents(i, :) = -cotangents(i, :);
    Laplacian_cotangents(i,i) = normalize_value2;
end



% solve x'Lx

% figure
% pdemesh(gm)
% hold on
% plot(gm.Mesh.Nodes(1,boundaryPts),gm.Mesh.Nodes(2,boundaryPts),"or","MarkerFaceColor","g")

a = zeros(numnodes, 1);
C = zeros(numnodes, numnodes);
for i =1:size(boundaryPts, 2)
    C(boundaryPts(i), boundaryPts(i)) = 1;
end

num_cagepts = size(cagePts, 2);
for i = 1:num_cagepts
    next = mod(i, num_cagepts)+1;
    C(cagePts(i), cagePts(i)) = -1;
    C(cagePts(i), cagePts(next)) = 1;
end

boundaryPts_loc = gm.Mesh.Nodes(:, boundaryPts);

%%% -------------------for multiple direction---------------------------- %%%
u_x = 6;
u_y = 6;
u_directions = [u_x u_y;
                u_x -u_y;
                -u_x -u_y;
                -u_x u_y];

for i=1:4
    boundaryPts_constraints = u_directions(i,:)*boundaryPts_loc;
    a(boundaryPts) = boundaryPts_constraints;
    if i > 1
        sol = max(sol, quadprog(Laplacian_cotangents,zeros(numnodes, 1),[],[],C,a));
    else
        sol = quadprog(Laplacian_cotangents,zeros(numnodes, 1),[],[],C,a);
    end
end
%%% -----------------------------------------------------------------%%%

%%% -------------------for one direction---------------------------- %%%
% u = [1 1];
% boundaryPts_constraints = u*boundaryPts_loc;
% a(boundaryPts) = boundaryPts_constraints;
% sol = quadprog(Laplacian_cotangents,zeros(numnodes, 1),[],[],C,a);
%%% -----------------------------------------------------------------%%%


patch_x = []; %holds the x values of the triangle vertices
patch_y = []; %holds the y values of the triangle vertices
patch_color = [];
face_gradient = [];
for i = 1:numelements
    v1_id = gm.Mesh.Elements(1, i);v2_id = gm.Mesh.Elements(2, i);v3_id = gm.Mesh.Elements(3, i);
    v1 = gm.Mesh.Nodes(:, v1_id);v2 = gm.Mesh.Nodes(:, v2_id);v3 = gm.Mesh.Nodes(:, v3_id);
    patch_x = [patch_x, [v1(1);v2(1);v3(1)]];
    patch_y = [patch_y, [v1(2);v2(2);v3(2)]];
    
    A = [(v1-v2); 0];
    B = [(v1-v3); 0];
    twoAT = norm(cross(A, B));
    
    B2 = (v1-v3)/twoAT;B2 = [-B2(2), B2(1)];
    B3 = (v2-v1)/twoAT;B3 = [-B3(2), B3(1)];

    grad_face = (sol(v2_id) - sol(v1_id))*B2 + (sol(v3_id)-sol(v1_id))*B3;

    face_gradient = [face_gradient; norm(grad_face)];
    patch_color = [patch_color, [sol(v1_id); sol(v2_id); sol(v3_id)]];
end

%%% -------------------setting up grid---------------------------- %%%
x = linspace(-3, 3, 400); y = linspace(-3, 3, 400);
[xGrid, yGrid] = meshgrid(x, y);
% https://www.mathworks.com/help/matlab/ref/scatteredinterpolant.html
F = scatteredInterpolant((gm.Mesh.Nodes(1,1:numnodes))', (gm.Mesh.Nodes(2,1:numnodes))', sol, 'linear');
tri_to_grid = F(xGrid, yGrid);
[grad_xx, grad_yy] = gradient(tri_to_grid, 6/100, 6/100);
magFX_grid = sqrt(grad_xx.^2 + grad_yy.^2);

exp_quantiles = quantile(magFX_grid, [0.025, 0.75], "all");
toobig = magFX_grid>exp_quantiles(2);
magFX_grid(toobig) = exp_quantiles(2);

bottom_10_percentile = quantile(magFX_grid, 0.08, "all");
bottom_10_percentile
inside_outside = double(magFX_grid <= bottom_10_percentile);

figure;
imagesc(magFX_grid);
%quiver(grad_xx, grad_yy);
axis xy; 
colorbar;

% n = 10;
% kernel = ones(n)/n.^2;
% inside_outside_conv = conv2(inside_outside, kernel,'same') ;

figure;
imagesc(inside_outside);
%quiver(grad_xx, grad_yy);
axis xy; 
colorbar;

figure;
contour(xGrid, yGrid, inside_outside, [0.5 0.5], 'k', 'LineWidth', 2);




%%% -----------------------------------------------------------------%%%


% plotting 
% patch for linear interpolation
%pdemesh(gm, "NodeLabels","on");
%axis equal

% imagesc([grad_xx, grad_yy]);
% pdemesh(gm);
exp_quantiles = quantile(face_gradient, [0.025, 0.75], "all");
toobig = face_gradient>exp_quantiles(2);
face_gradient(toobig) = exp_quantiles(2);

bottom_10_percentile = quantile(face_gradient, 0.1);
bottom_10_percentile
inside_outside = double(face_gradient <= bottom_10_percentile);


figure
patch(patch_x, patch_y, face_gradient)
% patch(patch_x, patch_y, patch_color)
colorbar

figure
patch(patch_x, patch_y, patch_color)
colorbar


figure
patch(patch_x, patch_y, inside_outside)
colorbar



vertex_values = zeros(numnodes, 1); % Initialize vertex values
vertex_counts = zeros(numnodes, 1); % To keep track of how many faces contribute to each vertex

for i = 1:size(inside_outside)
    face = gm.Mesh.Elements(1:3, i);
    vertex_values(face) = vertex_values(face) + inside_outside(i); % Accumulate face value for vertices
    vertex_counts(face) = vertex_counts(face) + 1; % Count contributions to each vertex
end
vertex_values = vertex_values ./ vertex_counts;

x_2 = linspace(-3, 3, 150); y_2 = linspace(-3, 3, 150);
[xGrid_2, yGrid_2] = meshgrid(x_2, y_2);
F_2 = scatteredInterpolant((gm.Mesh.Nodes(1,1:numnodes))', (gm.Mesh.Nodes(2,1:numnodes))', vertex_values, 'linear');
inside_grid = F(xGrid_2, yGrid_2);

figure;
contour(xGrid_2, yGrid_2, inside_grid, [0.5 0.5], 'k', 'LineWidth', 2);
colorbar;





