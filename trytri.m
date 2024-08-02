n = 5;r = 0.1;
c = exp(2i*pi*(1:n)/n);
rr = r*ones(size(c));
N = max(0, round(4+.5*log10(r)));
npts=3*N+2;
cagePts_num = n*npts;
xres = 5;
yres = 5;

x = linspace(-3, 3, 10); y = linspace(-3, 3, 10);


z = exp(pi*1i*(-50:50)'/50);
% x = [0.25, 0.75, 0.75, 0.25];
% y = [0.6, 0.6, 1.2, 1.2];
x = [-3, 3, 3, -3];
y = [-3, -3, 3, 3];
pgon = polyshape({x, real(c(1)+rr(1)*z)}, {y, imag(c(1)+rr(1)*z)});
for j = 2:n
    pgon = addboundary(pgon,real(c(j)+rr(j)*z),imag(c(j)+rr(j)*z));
end
tr = triangulation(pgon);
points = tr.Points; %row number is vertex id
connectivity = tr.ConnectivityList; %row number is triangle id, rows contain vertex id
% neighbors = neighbors(tr, id); %row is triangle id, tows contain triangle ids
gm = fegeometry(tr);
gm = generateMesh(gm, Hmin= 0.001);
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
  v1_id = gm.Mesh.Elements(1, i);v2_id = gm.Mesh.Elements(2, i);v3_id = gm.Mesh.Elements(3, i);

  v1 = gm.Mesh.Nodes(:, v1_id);v2 = gm.Mesh.Nodes(:, v2_id);v3 = gm.Mesh.Nodes(:, v3_id);
  
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

Laplacian_umbrella = zeros(numnodes);
Laplacian_cotangents = zeros(numnodes);
for i = 1:numnodes
    normalize_value1 = sum(neighborsMatrix(i, :));
    normalize_value2 = sum(cotangents(i, :));
    Laplacian_umbrella(i, :) = -neighborsMatrix(i, :)/normalize_value1;
    Laplacian_umbrella(i:i) = 1;
    Laplacian_cotangents(i, :) = -cotangents(i, :)/normalize_value2;
    Laplacian_cotangents(i,i) = 1;
end

Laplacian_cotangents


% solve x'Lx
cagePts
% figure
% pdemesh(gm)
% hold on
% plot(gm.Mesh.Nodes(1,boundaryPts),gm.Mesh.Nodes(2,boundaryPts),"or","MarkerFaceColor","g")

a = zeros(numnodes, 1);
C = zeros(numnodes, numnodes);
boundaryPts
for i =1:size(boundaryPts, 2)
    C(boundaryPts(i), boundaryPts(i)) = 1;
end

num_cagepts = size(cagePts, 2);
for i = 1:num_cagepts
    next = mod(i, num_cagepts)+1;
    C(cagePts(i), cagePts(i)) = -1;
    C(cagePts(i), cagePts(next)) = 1;
end
C

boundaryPts_loc = gm.Mesh.Nodes(:, boundaryPts);
u = [1 1];
boundaryPts_constraints = u*boundaryPts_loc;

%a_boundary = a(boundaryPts);

a(boundaryPts) = boundaryPts_constraints;

x = quadprog(Laplacian_cotangents,zeros(numnodes, 1),[],[],C,a);
x
patch_x = []
patch_y = []
patch_color = []
for i = 1:numelements
    v1_id = gm.Mesh.Elements(1, i);v2_id = gm.Mesh.Elements(2, i);v3_id = gm.Mesh.Elements(3, i);
    v1 = gm.Mesh.Nodes(:, v1_id);v2 = gm.Mesh.Nodes(:, v2_id);v3 = gm.Mesh.Nodes(:, v3_id);
    patch_x = [patch_x, [v1(1);v2(1);v3(1)]];
    patch_y = [patch_y, [v1(2);v2(2);v3(2)]];
    patch_color = [patch_color, [abs(x(v1_id)); abs(x(v2_id)); abs(x(v3_id))]];
end
% plotting 
% patch for linear interpolation
%pdemesh(gm, "NodeLabels","on");
%axis equal
figure
patch(patch_x, patch_y, patch_color)
colorbar
