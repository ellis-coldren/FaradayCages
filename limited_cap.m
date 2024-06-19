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
disp(centers);

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
disp(size(P));
disp(size(rhs));
X = P\rhs;
disp(X(1));
e = X(1); X(1) = []; %constant voltage on wires
d = X(1:2*N+1:end); X(1:2*N+1:end) = [];
disp(d(1));
a = X(1:2:end); b=X(2:2:end);
disp(size(X));

x = linspace(-1.4, 2.2, 120); y = linspace(-1.8, 1.8, 120);
[xx, yy] = meshgrid(x, y); zz=xx+1i*yy; uu=log(abs(zz-zs));

V = zeros(size(zz));

% Calculate potential at each point in the grid
for i = 1:numel(zz)
    point = zz(i);
    for j = 1:N
        ri = centers(j);  % Position of point charge j
        distance = norm(point - ri);  % Distance between r and ri
        V(i) = V(i) + q_lc(j) / distance;  % Sum contributions from all point charges(wires)
    end
    V(i) = V(i)-log(abs(point-zs)); %add contribution from the external field
end
disp(centers);

figure;
contourf(xx, yy, abs(V), 50);
hold on;
%quiver(real(V), imag(V));
for j=1:n, disk = centers(j)+rr(j)*z; fill(real(disk), imag(disk), [1 .7 .7])
    hold on, plot(disk, '-r'), end
colorbar;
colormap jet; 