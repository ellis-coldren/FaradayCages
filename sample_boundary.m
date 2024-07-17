n = 20;r = 0.1;
c = exp(2i*pi*(1:n)/n);
rr = r*ones(size(c));
N = max(0, round(4+.5*log10(r)));
npts=3*N+2;
npts=3;

disp(size(c));
%location of the boundary
global min_x max_x min_y max_y;
min_x = -1.5;
max_x = 1.5;
min_y = -1.5;
max_y = 1.5;
dist_x = max_x - min_x;
dist_y = max_y - min_y;
setup_x_boundary = min_x:dist_x/(n-1):max_x;
setup_y_boundary = min_y:dist_y/(n-1):max_y;
top_boundary = setup_x_boundary' + (max_y-r)*1i;
bottom_boundary = setup_x_boundary' + (min_y+r)*1i;
right_boundary = (min_x+r) + setup_y_boundary'*1i;
left_boundary = (max_x-r) + setup_y_boundary'*1i;

boundary_sample = [top_boundary; bottom_boundary; right_boundary; left_boundary];

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

A = [0; -ones(size(z, 1), 1); zeros(4*size(top_boundary, 1), 1)];
zs = 2;
%explicitly define the values on boundary
linear_vector = [4 6];
top_boundary_values = zeros(size(top_boundary, 1), 1);
bottom_boundary_values = zeros(size(top_boundary, 1), 1);
right_boundary_values = zeros(size(top_boundary, 1), 1);
left_boundary_values = zeros(size(top_boundary, 1), 1);
for j = 1:size(top_boundary, 1)
    top_boundary_values(j) = linear_vector * [real(top_boundary(j)); imag(top_boundary(j))];
    bottom_boundary_values(j) = linear_vector * [real(bottom_boundary(j)); imag(bottom_boundary(j))];
    right_boundary_values(j) = linear_vector * [real(right_boundary(j)); imag(right_boundary(j))];
    left_boundary_values(j) = linear_vector * [real(left_boundary(j)); imag(left_boundary(j))];
end

boundary_values = [top_boundary_values; bottom_boundary_values; right_boundary_values; left_boundary_values];

rhs = [0; zeros(size(z, 1), 1); top_boundary_values; bottom_boundary_values; right_boundary_values; left_boundary_values];
% rhs = [0; -log(abs(z-zs)); top_boundary_values];

disp(top_boundary');

disp(abs(top_boundary(2) - top_boundary(1)));
for j=1:n
    A = [A [1; log(abs(z-c(j))); log(abs(top_boundary-top_boundary(j))); log(abs(bottom_boundary-bottom_boundary(j))); log(abs(right_boundary - right_boundary(j))); log(abs(left_boundary - left_boundary(j)))]];
    A(A== -Inf) = 0;
    for k = 1:N
        zck = (z-c(j)).^(-k);
        disp(((top_boundary-top_boundary(j)).^(-k))');
        zck_top = (top_boundary-top_boundary(j)).^(-k);
        zck_top(zck_top== Inf) = 0;
        zck_bottom = (bottom_boundary-bottom_boundary(j)).^(-k);
        zck_bottom(zck_bottom== Inf) = 0;
        zck_right = (right_boundary - right_boundary(j)).^(-k);
        zck_right(zck_right== Inf)=0;
        zck_left = (left_boundary - left_boundary(j)).^(-k);
        zck_left(zck_left== Inf) = 0;
        A = [A [0; real(zck); real(zck_top); real(zck_bottom); real(zck_right); real(zck_left)] [0; imag(zck); imag(zck_top); imag(zck_bottom); imag(zck_right); imag(zck_left)]];
    end
end

disp(A);
X = A\rhs;
e=X(1); X(1) = []; %constant voltage on wires
disp(X);
d = X(1:2*N+1:end); X(1:2*N+1:end) = [];
a = X(1:2:end); b=X(2:2:end);


x = linspace(-1.4, 2.2, 120); y = linspace(-1.8, 1.8, 120);
[xx, yy] = meshgrid(x, y); zz=xx+1i*yy; uu=zeros(size(zz));
% disp(uu)


for j=1:n
    uu = uu+d(j)*log(abs(zz-c(j)) + abs(zz-top_boundary(j)) + abs(zz-bottom_boundary(j)) + abs(zz-right_boundary(j)) + abs(zz-left_boundary(j)));
    for k=1:N, zck = (zz-c(j) + zz-top_boundary(j) + zz-bottom_boundary(j) + zz-right_boundary(j)+zz-left_boundary(j)).^(-k); kk=k+(j-1)*N;
        uu = uu+a(kk)*real(zck)+b(kk)*imag(zck); end
end

for j=1:n, uu(abs(zz-c(j))<rr(j)) = NaN; end
z = exp(pi*1i*(-50:50)'/50);

[grad_xx, grad_yy] = gradient(real(uu), 3.6/120, 3.6/120);
magFX_grid = sqrt(grad_xx.^2 + grad_yy.^2);
exp_quantiles = quantile(magFX_grid, [0.025, 0.75], "all");
disp(exp_quantiles);
toobig = magFX_grid>exp_quantiles(2);
magFX_grid(toobig) = exp_quantiles(2);


for j=1:n, disk = c(j)+rr(j)*z; fill(real(disk), imag(disk), [1 .7 .7])
    hold on, plot(disk, '-r'), end
contour(xx, yy, real(uu), -2:.1:1.2), colormap([0 0 0]), axis([-1.4 2.2 -1.8 1.8])
axis square, plot(real(zs), imag(zs), '.r')
figure;
    %quiver(grad_xx, grad_yy);
    imagesc(magFX_grid);
colorbar;
axis equal;




