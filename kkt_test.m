
n = 2;r = 0.01; % frequency and radius of cage vertices
t = linspace(0, 2*pi, 1000); % parameter value
using_vec = 1;

% -------------------- Circle ----------------------- %
radius = 0.75;
x_shape = radius*cos(t);
y_shape = radius*sin(t);

% x_shape = [x_shape radius*cos(t) + 0.5];
% y_shape = [y_shape radius*sin(t)];

% -------------------- Heart ----------------------- %
% x_shape = 0.08*(16*sin(t).^3);
% y_shape = 0.08*(13*cos(t) - 5*cos(2*t) - 2*cos(3*t) - cos(4*t));

% -------------------- Diamond ----------------------- %
% a = 0.7;
% x_shape = a*cos(t).^3;
% y_shape = a*sin(t).^3;

% x_shape = [x_shape a*cos(t).^3 + 0.5];
% y_shape = [y_shape a*sin(t).^3];

% -------------------- Rose curve ----------------------- %
% a = 0.6;
% petals = 4;
% x_shape = a*cos(petals.*t).*cos(t);
% y_shape = a*cos(petals.*t).*sin(t);

% -------------------- Hypotrochoid ----------------------- %
%(sometimes playing around messes up, see https://www.desmos.com/calculator/r8bsgvsszo)

% fr = 1.6;
% rr = 0.2;
% d = 0.5;
% x_shape = 0.5*((fr-rr)*cos(t)+d*cos(((fr-rr)/rr)*t));
% y_shape = 0.5*((fr-rr)*sin(t) - d*sin(((fr-rr)/rr)*t));

% -------------------- Reading .vec files ----------------------- %
% vecFile = 'fruit.vec';
% vecRead = xmlread(vecFile);
% using_vec = 0;
% 
% sketch_edges = vecRead.getElementsByTagName('edge');
% x_shape = [];
% y_shape = [];
% 
% for i = 0:sketch_edges.getLength-1
%     edge = sketch_edges.item(i);
%     % the xywdense has format: xywdense(layer(int) x1, y1, w1, x2, y2, w2, ...);
%     densesample = char(edge.getAttribute('curve'));
%     sample_substring = regexp(densesample, '[+-]?\d*\.?\d+', 'match');
%     sample_double = str2double(sample_substring);
% 
%     length_x = length(sample_double(2:3:end));
%     length_y = length(sample_double(3:3:end));
% 
%     x_indices = round(linspace(1, length_x, n));
%     y_indices = round(linspace(1, length_y, n));
%     
%     x_sample = sample_double(2:3:end);
%     x_to_add = x_sample(x_indices);
%     y_sample = sample_double(3:3:end);
%     y_to_add = y_sample(y_indices);
% 
%     x_shape = [x_shape, x_sample];
%     y_shape = [y_shape, y_sample];
% end
% min_x = min(x_shape);
% max_x = max(x_shape);
% min_y = min(y_shape);
% max_y = max(y_shape);
% x_scale = (3 - (-3))/(max_x - min_x);
% y_scale = (3 - (-3))/(max_y - min_y);
% x_shape = x_scale*x_shape;
% y_shape = -y_scale*y_shape;

% --------------------------------------------------------------- %



noise_param = 0;  % param to control deg of randomness in curve

% TODO normrnd - matlab gaussian model for sampling
p = [x_shape', y_shape'];
q = curvspace(p,n+1); % generates points that interpolate curve
xx = q(:, 1) + normrnd(0, 1, [size(q, 1), 1])*noise_param; % x values of curve
yy = q(:, 2) + normrnd(0, 1, [size(q, 1), 1])*noise_param; % y values of curve
if using_vec == 0
    xx = p(:, 1)+ normrnd(0, 1, [size(p, 1), 1])*noise_param;
    yy = p(:, 2) + normrnd(0, 1, [size(p, 1), 1])*noise_param;;
end
c = xx + 1i*yy; % curve as vector of complex coordinates
rr = r*ones(size(c)); % radii of cage vertices as vector
N = max(0, round(4+.5*log10(r)));
npts=3*N+2;
% cagePts_num = n*npts;
% xres = 5;
% yres = 5;

z = exp(pi*1i*(-10:10)'/10);
% z = 1;
% z = [1 -1]; % vector of circle coordinates for cage vertices
min_x = min(x_shape);
max_x = max(x_shape);
min_y = min(y_shape);
max_y = max(y_shape);
% x = [min_x-2, max_x+2, max_x+2, min_x-2]; % boundary of region
% y = [min_y-2, min_y-2, max_y+2, max_y+2]; % boundary of region

x = [-1, 1, 1, -1]; % boundary of region
y = [-1, -1, 1, 1]; % boundary of region

% https://www.mathworks.com/help/matlab/ref/polyshape.html %
% ---------------- setting up cage vertices --------------------- %

% connects vertices of x and y, adds hole where first cage vertex is (c(1) is the center of the cage point, + rr(1)*z draws a circle with radius rr(1) around it) %
pgon = polyshape({x, real(c(1)+rr(1)*z)}, {y, imag(c(1)+rr(1)*z)}); 
for j = 2:length(c)-1
    % adds another hole for each cage vertex %
    pgon = addboundary(pgon,real(c(j)+rr(j)*z),imag(c(j)+rr(j)*z));
end

% pgon = polyshape({x, real(c(1))}, {y, imag(c(1))}); 
% for j = 2:length(c)-1
%     % adds another hole for each cage vertex %
%     pgon = addboundary(pgon,real(c(j)),imag(c(j)));
% end

% triangulates the region %
tr = triangulation(pgon);
points = tr.Points; %row number is vertex id
connectivity = tr.ConnectivityList; %row number is triangle id, rows contain vertex id

% creates finite element model so that we can create a mesh %
% https://www.mathworks.com/help/pde/ug/fegeometry.html

gm = fegeometry(tr);
gm = generateMesh(gm, Hmax = 2, Hmin= 0.01);
figure;
pdemesh(gm, 'NodeLabels','on');
numelements = size(gm.Mesh.Elements, 2);
numnodes = size(gm.Mesh.Nodes, 2);
        % elements is 6 x (number of triangles)
        % the vertex ids of the triangles are in the first three rows of gm.Mesh.Elements
        % midpoints of the edges of the triangle are the last 3 rows of the gm.Mesh.Elements
        % vertices of the element(triangle) are the first 3 rows of gm.Mesh.Elements

% this finds which vertices of our mesh lie on a cage vertex:
cagePts = findNodes(gm.Mesh,"radius",[real(c(1)), imag(c(1))], rr(1));
for j = 1:length(c)
    cagePts = [cagePts findNodes(gm.Mesh,"radius",[real(c(j)), imag(c(j))], rr(j))];
end

% this finds which vertices of our mesh lie on the boundary:
boundaryPts = findNodes(gm.Mesh, "region", "Edge",[1 2 3 4]);

neighborsMatrix = sparse(numnodes, numnodes);
%neighborsMatrix(i, j) = neighborsMatrix(j, i) == 1 if vertex i and j share an edge, and 0 otherwise
cotangents = sparse(numnodes, numnodes);
% cotangents(i, j) = 1/2*cot(alpha(i, j)) + cot(beta(i, j))
for i = 1:numelements %iterate through the triangles
  v1_id = gm.Mesh.Elements(1, i); v2_id = gm.Mesh.Elements(2, i); v3_id = gm.Mesh.Elements(3, i); %get vertex id's of the triangle
  fprintf('triangle = %d %d %d \n', v1_id, v2_id, v3_id);
  v1 = gm.Mesh.Nodes(:, v1_id);v2 = gm.Mesh.Nodes(:, v2_id);v3 = gm.Mesh.Nodes(:, v3_id); % get vertex locations
  
  neighborsMatrix(v1_id, v2_id) = 1;neighborsMatrix(v2_id, v3_id) = 1;neighborsMatrix(v3_id, v1_id) = 1;
  neighborsMatrix(v2_id, v1_id) = 1;neighborsMatrix(v3_id, v2_id) = 1;neighborsMatrix(v1_id, v3_id) = 1;


  cotangents(v1_id, v2_id) = cotangents(v1_id, v2_id) + 0.5*cot(acos(dot(v1-v3, v2-v3)/(norm(v1-v3)*norm(v2-v3))));
  cotangents(v2_id, v1_id) = cotangents(v1_id, v2_id);
  fprintf('v1 v2 cotangent \n', full(cotangents(v1_id, v2_id)));
  
  cotangents(v2_id, v3_id) = cotangents(v2_id, v3_id) + 0.5*cot(acos(dot(v2-v1, v3-v1)/(norm(v2-v1)*norm(v3-v1))));
  cotangents(v3_id, v2_id) = cotangents(v2_id, v3_id);
  fprintf('v3 v2 cotangent \n', full(cotangents(v2_id, v3_id)));
  
  cotangents(v3_id, v1_id) = cotangents(v3_id, v1_id) + 0.5*cot(acos(dot(v3-v2, v1-v2)/(norm(v3-v2)*norm(v1-v2))));
  cotangents(v1_id, v3_id) = cotangents(v3_id, v1_id);
  fprintf('v1 v3 cotangent \n', full(cotangents(v3_id, v1_id)));

end

% if row or column i is entirely 0, then vertex i is a midpoint
%%%% NOTE: if vertex i is a midpoint, than any vertex with an id greater
%%%% than i is also a midpoint %%%%
% the following removes those midpoints for efficiency
neighborsMatrix = neighborsMatrix(any(neighborsMatrix, 2), :); %removing rows of the middle point nodes
neighborsMatrix = neighborsMatrix(:, any(neighborsMatrix, 1)); %removing columns of the middle point nodes
cotangents = cotangents(any(cotangents, 2), :); %removing rows of the middle point nodes
cotangents = cotangents(:, any(cotangents, 1)); %removing columns of the middle point nodes

%updates numnodes to exclude midpoints
numnodes = size(neighborsMatrix, 1);

% updates cagePts and boundaryPts to remove any midpoint vertices
colstodelete = cagePts > numnodes;
cagePts(colstodelete) = [];
colstodelete = boundaryPts > numnodes;
boundaryPts(colstodelete) = [];

% creating cotangent Laplacian
Laplacian_cotangents = sparse(numnodes, numnodes);
for i = 1:numnodes
    normalize_value2 = sum(cotangents(i, :));
    Laplacian_cotangents(i, :) = -cotangents(i, :);
    Laplacian_cotangents(i,i) = normalize_value2;
end


C = sparse(numnodes, numnodes);
% C(i, i) = 1 if vertex i is on the boundary
for i =1:size(boundaryPts, 2)
    C(boundaryPts(i), boundaryPts(i)) = 1;
end

num_cagepts = size(cagePts, 2);
% if i is a cage vertex
% get the next cage vertex
% C(i, i) = -1 and C(i, next) = 1
% so that all cage vertices are equal
for i = 1:(num_cagepts-1)
    next = mod(i, num_cagepts)+1;
    C(cagePts(i), cagePts(i)) = -1;
    C(cagePts(i), cagePts(next)) = 1;
end
C = C(any(C, 2), :);
a = zeros(size(C, 1), 1);

% coordinates of the boundary points
boundaryPts_loc = gm.Mesh.Nodes(:, boundaryPts);

%%% -------------------for multiple direction---------------------------- %%%
% num_direcs = 15;    % num directions to test
% fieldvec_mag = 1;
% u_directions = [];
% for j = 1:num_direcs
%     u_directions = [u_directions; fieldvec_mag*cos(j*2*pi/num_direcs) fieldvec_mag*sin(j*2*pi/num_direcs)];
% end
% sol = zeros(numnodes, num_direcs);
% for i=1:num_direcs
%     boundaryPts_constraints = u_directions(i,:)*boundaryPts_loc;
%     a(boundaryPts) = boundaryPts_constraints;
%     sol(:, i) = quadprog(Laplacian_cotangents,zeros(numnodes, 1),[],[],C,a);
% end
% size(sol)
%%% -----------------------------------------------------------------%%%

%%% -------------------for one direction---------------------------- %%%
u = [1 1];
num_direcs = 1;
boundaryPts_constraints = u*boundaryPts_loc;
a(boundaryPts) = boundaryPts_constraints;
sol = quadprog(Laplacian_cotangents,zeros(numnodes, 1),[],[],C,a);


%%% -------------------different KKT Solvers---------------------------- %%%
LHS = [Laplacian_cotangents, C.'; 
      C, sparse(size(C, 1), size(C, 1))];
RHS = [zeros(numnodes, 1); a];
x_v = LHS\RHS;
sol_2 = x_v(1:numnodes);

x_v = lsqr(LHS, RHS, 1e-6, 2000);
sol_3 = x_v(1:numnodes);

x_v = pcg(LHS, RHS);
sol_4 = x_v(1:numnodes);


a_2 = zeros(numnodes, 1);
C_2 = sparse(numnodes, numnodes);
num_cagepts = size(cagePts, 2);
for i = 1:num_cagepts
    next = mod(i, num_cagepts)+1;
    C_2(cagePts(i), cagePts(i)) = -1;
    C_2(cagePts(i), cagePts(next)) = 1;
end
sol_5 = quadprog(Laplacian_cotangents,zeros(numnodes, 1),[],[],C_2,a_2);

%%% -----------------------------------------------------------------%%%


patch_x = []; %holds the x values of the triangle vertices
patch_y = []; %holds the y values of the triangle vertices
patch_color = []; %holds the interpolated value of sol of the patch
patch_color_2 = [];
patch_color_3 = [];
patch_color_4 = [];
patch_color_5 = [];
face_gradient = []; % holds the magnitude of the face gradient of the patch
%disp('LHS = '); disp(full(LHS));
disp('Cage = '); disp(size(C));
%disp('Laplacian_cot = '); disp(full(Laplacian_cotangents));
[V,val] = eigs(LHS,1,'smallestabs');
disp('eigenvector ='); disp(V);
disp('eigenvalue = '); disp(val);
disp(size(LHS));
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

    for i = 1:num_direcs
        grad_face = (sol(v2_id, i) - sol(v1_id, i))*B2 + (sol(v3_id, i)-sol(v1_id, i))*B3;
        if i == 1
            norm_grad = norm(grad_face);
        else
            norm_grad = max(norm_grad, norm(grad_face));
        end
    end

    face_gradient = [face_gradient; norm_grad];
    patch_color = [patch_color, [sol(v1_id, 1); sol(v2_id, 1); sol(v3_id, 1)]];
    patch_color_2 = [patch_color_2, [sol_2(v1_id, 1); sol_2(v2_id, 1); sol_2(v3_id, 1)]];
    patch_color_3 = [patch_color_3, [sol_3(v1_id, 1); sol_3(v2_id, 1); sol_3(v3_id, 1)]];
    patch_color_4 = [patch_color_4, [sol_4(v1_id, 1); sol_4(v2_id, 1); sol_4(v3_id, 1)]];
    patch_color_5 = [patch_color_5, [sol_5(v1_id, 1); sol_5(v2_id, 1); sol_5(v3_id, 1)]];
end

% disp([gm.Mesh.Nodes(:, boundaryPts)', sol(boundaryPts)]);
exp_quantiles = quantile(face_gradient, [0.025, 0.75], "all");
toobig = face_gradient>exp_quantiles(2);
face_gradient(toobig) = exp_quantiles(2);

bottom_15_percentile = quantile(face_gradient, 0.25);
inside_outside = double(face_gradient <= bottom_15_percentile);

max(patch_color, [], 'all');
figure;
patch(patch_x, patch_y, patch_color);
title('QuadProg Solve');
colorbar;

figure;
patch(patch_x, patch_y, patch_color_2);
title('KKT Solve (LHS\RHS)');
colorbar;

figure;
patch(patch_x, patch_y, patch_color_3);
title('KKT Solve (lsqr)');
colorbar;

figure;
patch(patch_x, patch_y, patch_color_4);
title('KKT Solve (pcg)');
colorbar;

figure;
patch(patch_x, patch_y, patch_color_5);
title('No Boundary');
colorbar;


% figure;
% patch(patch_x, patch_y, face_gradient);
% title('Mag Grad');
% colorbar;
% 
% 
% figure;
% patch(patch_x, patch_y, inside_outside);
% title('Inside/Outside');
% colorbar;
% 
% 
% 
% vertex_values = zeros(numnodes, 1); % Initialize vertex values
% vertex_counts = zeros(numnodes, 1); % To keep track of how many faces contribute to each vertex
% 
% for i = 1:size(inside_outside)
%     face = gm.Mesh.Elements(1:3, i);
%     vertex_values(face) = vertex_values(face) + face_gradient(i); % Accumulate face value for vertices
%     vertex_counts(face) = vertex_counts(face) + 1; % Count contributions to each vertex
% end
% vertex_values = vertex_values ./ vertex_counts;
% 
% figure;
% tricontf(gm.Mesh.Nodes(1,1:numnodes)', gm.Mesh.Nodes(2,1:numnodes)', gm.Mesh.Elements(1:3, :)', vertex_values, [0.5 0.5], 'k');





