% linear potential boundary conditions

% Step 1:
% setup electric potential that varies linearly along the boundary
a = 0.5;
b = 1.2;
%for x, plug in points on the cage
% for x from [-2, 2]
% for y from [-2i, 2i]
% setting up portion of grid
x = linspace(-1.4, 2.2, 120); 
y = linspace(-1.8, 1.8, 120);

uu = zeros(120, 120); %setting up potential field
%left boundary, boundary condition is linear
uu(:, 1) = y';
%right boundary
uu(:, end) = y' + (2.2+1.4);
%top boundary
uu(1, :) = x;
%bottom boundary
uu(end,:) = x+(1.8+1.8);
%solve Laplace's equation?
% https://www.public.asu.edu/~hhuang38/pde_slides_numerical_laplace.pdf
% Gauss-Seidel for matrix problem

tolerance = 1e-5;
max_iter = 1000;

for iter = 1:max_iter
    uu_old = uu;
    for i=2:119
        for j=2:119
            uu(i, j) = 0.25 * (uu(i+1, j) + uu(i, j+1) + uu(i, j+1) + uu(i, j-1));
        end
    end

    if max(max(abs(uu - uu_old))) < tolerance
        break;
    end
end

[grad_xx, grad_yy] = gradient(real(uu), 3.6/120, 3.6/120);
[xx, yy] = meshgrid(x, y); 
contour(xx, yy, real(uu), -2:.1:1.2), colormap([0 0 0]), axis([-1.4 2.2 -1.8 1.8]);
figure;
    quiver(grad_xx, grad_yy);
    %imagesc(magFX_grid);
colorbar;
axis equal

