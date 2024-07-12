% linear potential boundary conditions

% Step 1:
% setup electric potential that varies linearly along the boundary
a = 0.5;
b = 1.2;
%for x, plug in points on the cage
% for x from [-2, 2]
% for y from [-2i, 2i]
% setting up portion of grid
N = 120; % number of steps
h = 3/N; % step size
x = -1.5:h:1.5; 
y = -1.5:h:1.5;
[X, Y] = meshgrid(x, y);

w = zeros(N+1, N+1); 
boundary_vec = [1, 1];
% setting up boundary:
function a = myboundary(x,y, boundary_vec)
    a = boundary_vec*[x, y]';
end
for k=1:N %vertical boundaries
    leftboundary = myboundary((k-1)*h, 0, boundary_vec);
    rightboundary = myboundary((k-1)*h, (N-1), boundary_vec);
    w(k, 1) = leftboundary;
    w(k, N) = rightboundary;
end

for j=0:N %horizontal boundaries
    lowerboundary = myboundary(0, (k-1)*h, boundary_vec);
    upperboundary = myboundary(N-1, (k-1)*h, boundary_vec);
    w(1, k) = leftboundary;
    w(N, k) = rightboundary;
end

%now, write as system of (N-1)x(N-1)equations Aw = r
N2 = (N-1)*(N-1);
A = zeros(N2, N2);
disp(size(A));
%setup diagonal
for i=1:N-1
    for j=1:N-1
        A((i-1)+(N-1)*(j-1)+1, (i-1)+(N-1)*(j-1)+1)=-4;
    end
end
%lower diag, should be ones
for i=2:N-1
    for j=1:N-1
        A((i-1)+(N-1)*(j-1)+1, (i-1)+(N-1)*(j-1)) = 1;
    end
end

%upper diag
for i=1:N-2
    for j=1:N-1
        A((i-1)+(N-1)*(j-1)+1, (i-1)+(N-1)*(j-1)+2) = 1;
    end
end

%lower id matrix
for i=1:N-1
    for j=2:N-1
        A((i-1)+(N-1)*(j-1)+1, (i-1)+(N-1)*(j-2)+1)=1;
    end
end

%lower id matrix
for i=1:N-1
    for j=1:N-2
        A((i-1)+(N-1)*(j-1)+1, (i-1)+(N-1)*(j)+1)=1;
    end
end

disp(size(A));
boundary_condit = zeros(1, N2);
for i=1:N-1
    subvector = zeros(1, N-1); %ith column
    for j=2:N-1
        subvector(j)=w(j, i);
    end
    boundary_condit((i-1)*(N-1)+1: i*(N-1)) = subvector';
end

% making r
bx = zeros(1, N2);
for j=1:N-1
    bx(j*(N-1)) = w(1, j);
    bx(j*(N-1)+N-2) = w(N, j);
end

by = zeros(1, N2);
subvector_1 = zeros(1, N-1);
subvector_2 = zeros(1, N-1);
for j=1:N-1
    subvector_1(j) = w(j, 1); 
    subvector_2(j) = w(j, N);
end
by(1:(N-1)) = subvector_1;
by(N2-(N-1)+1: N2) = subvector_2;

r = zeros(1, N2);
disp(N2);
disp(size(r));
for j=1:N-1
    r((j-1)*(N-1)+1:j*(N-1)) = -bx((j-1)*(N-1)+1:j*(N-1)) - by((j-1)*(N-1)+1:j*(N-1));
end

A_inv = inv(A);
disp(size(A_inv));
disp(size(r'));
C = A_inv*r';

boundary_condit(1:N-1, 1:N-1) = reshape(C, N-1, N-1);

[grad_xx, grad_yy] = gradient(real(w), 3/120, 3/120);
contour(X, Y, w, -2:.1:1.2), colormap([0 0 0]), axis([-1.5 1.5 -1.5 1.5]);
figure;
    quiver(grad_xx, grad_yy);
    %imagesc(magFX_grid);
colorbar;
axis equal

