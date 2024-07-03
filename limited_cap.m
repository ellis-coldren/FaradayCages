%needs optimization toolbox install

n = 80; %number of wires
r = 0.02; %radius of wires
t = linspace(0, 2*pi, 1000);

%finding centers of the disks
%Circle
radius = 0.75;
x = radius*cos(t);
y=radius*sin(t);

p = [x', y'];
samplepoints = curvspace(p,n);
xx = samplepoints(:, 1);
yy = samplepoints(:, 2);
centers = xx + 1i*yy;

zs = [1.5]; %location of external field

%setting up A
quad_vec = zeros(size(centers));
quad_vec(:) = -log(r);
A = diag(quad_vec);
for k=1:size(centers)
    quad_vec = zeros(size(centers));
    for j=(k+1):size(centers)
        quad_vec(j) = -log(abs(centers(k) - centers(j)));
    end
    quad_vec(k) = -log(r);
    A(k,:)=quad_vec;
end

%setting up f
f = zeros(size(centers));
for k=1:size(centers)
    f(k) = -log(abs(centers(k) - zs));
end

%setting up c
c = ones(size(centers));
x0 = ones(size(centers));
A_eq = c';
b_eq = 0;

%can't use quadprog, problem is not convex
% [q, fval] = quadprog(A, f, [], [], A_eq, b_eq);


func = @(x)quadconstreq(x, A, f);
x = fmincon(func, x0, [], [], A_eq, b_eq);
%fmincon :: problem appears unbounded
%fmincon stopped because the objective function value is less than the
%value of the obkective function limit, and constraints are satisfied to
%within the value of the constraint tolerance

%setting up quadratic solve
M = [A, c; c', 0];
rhs = [f; 0];

solution = M \ rhs;

q_lc = solution(1:n);
lambda = solution(n);
disp(q_lc);

% syms q [size(centers) 1]
% A_sym = sym(A); %have to use symbolic
% f_sym = sym(f);
% c_sym = sym(c);

% E(1) = (1/2)*q'*A*q-f'*q==0;
% E(2) = c'*q==0;
% 
% q_zc = solve(E, q, 'returnconditions', true);

rr = r*ones(size(centers));
N = max(0, round(4+.5*log10(r)));

npts = 3*N+2;
circ = exp((1:npts)'*2i*pi/npts);
z = []; for j=1:n
    z=[z;centers(j)+rr(j)*circ]; end
P = [0; -ones(size(z))];

%% allows for multiple point charges
rhs_log = 0;
for i = 1:width(zs)
    rhs_log = rhs_log - log(abs(z-zs(i)));
end
rhs = [0; -log(abs(z-zs))];
for j=1:n
    P = [P [1; log(abs(z-centers(j)))]];
    for k = 1:N
        zck = (z-centers(j)).^(-k);
        P = [P [0; real(zck)] [0;imag(zck)]];
    end
end

X = P\rhs;
e = X(1); X(1) = []; %constant voltage on wires
d = X(1:2*N+1:end); X(1:2*N+1:end) = [];
a = X(1:2:end); b=X(2:2:end);

x = linspace(-1.4, 2.2, 120); y = linspace(-1.8, 1.8, 120);
[xx, yy] = meshgrid(x, y); zz=xx+1i*yy; uu=log(abs(zz-zs));

V = zeros(size(zz));

% Calculate potential at each point in the grid
for i = 1:numel(zz)
    point = zz(i);
    for j = 1:N
        ri = centers(j);  % Position of point charge j
        %distance = norm(point - ri);  % Distance between r and ri
        V(i) = V(i) + log(abs(q_lc(j)-point));  % Sum contributions from all point charges(wires)
    end
    V(i) = V(i)+log(abs(point-zs)); %add contribution from the external field
end


% figure;
% quiver(xx, yy, V, 50);
% hold on;
% %quiver(real(V), imag(V));
% for j=1:n, disk = centers(j)+rr(j)*z; fill(real(disk), imag(disk), [1 .7 .7])
%     hold on, plot(disk, '-r'), end
% colorbar;
% colormap jet;

[grad_xx, grad_yy] = gradient(real(V), 3.6/120, 3.6/120);

%grad_xx(in)=0;grad_yy(in)=0; %filter for too close to point charge
magFX_grid = sqrt(grad_xx.^2 + grad_yy.^2);
%magFX_grid(in)=max(magFX_grid, [], "all"); %filter for too close to point charge

exp_quantiles = quantile(magFX_grid, [0.025, 0.65], "all");
toobig = magFX_grid>exp_quantiles(2);
magFX_grid(toobig) = exp_quantiles(2);

for j=1:n, V(abs(zz-centers(j))<rr(j)) = NaN; end
z = exp(pi*1i*(-50:50)'/50);
for j=1:n, disk = centers(j)+rr(j)*z; fill(real(disk), imag(disk), [1 .7 .7])
    hold on, plot(disk, '-r'), end
contour(xx, yy, real(V), -2:.1:2), colormap([0 0 0]), axis([-1.4 2.2 -1.8 1.8])
axis square, plot(real(zs), imag(zs), '.r')
figure;
    % quiver(grad_xx, grad_yy);
    % imagesc(magFX_grid);
colorbar;
axis equal
