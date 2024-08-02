n = 5;r = 0.1;
c = exp(2i*pi*(1:n)/n);
rr = r*ones(size(c));
N = max(0, round(4+.5*log10(r)));
npts=3*N+2;
cagePts = n*npts;
xres = 20;
yres = 20;

%location of the boundary
global min_x max_x min_y max_y;
min_x = -3;
max_x = 3;
min_y = -3;
max_y = 3;
dist_x = max_x - min_x;
dist_y = max_y - min_y;
setup_x_boundary = min_x:dist_x/(xres-1):max_x;
setup_y_boundary = min_y:dist_y/(yres-1):max_y;
top_boundary = setup_x_boundary' + (max_y)*1i;
bottom_boundary = setup_x_boundary' + (min_y)*1i;
right_boundary = (min_x) + setup_y_boundary'*1i;
left_boundary = (max_x) + setup_y_boundary'*1i;

boundary_sample = [top_boundary; bottom_boundary; right_boundary; left_boundary];
disp(boundary_sample);

% Not used
function d = dist_to_top(zx)
      global max_x min_x max_y
      pt = [real(zx), imag(zx), 0];
      a = [max_x, max_y, 0] - [min_x, max_y, 0];
      b = pt - [min_x, max_y, 0];
      d = norm(cross(a,b)) / norm(a);
end

circ=exp((1:npts)'*2i*pi/npts);
z = []; for j=1:n
    z=[z;c(j)+rr(j)*circ]; end
z = [z; boundary_sample];

P = [imag(z) real(z)];

x = linspace(-3, 3, 120); y = linspace(-3, 3, 120);

DT = delaunay(P);
triplot(DT, x, y);
