n = 100;
r = 0.02;
t = linspace(0, 2*pi, 1000);

%Circle
radius = 100000;
x = radius*cos(t);
y=radius*sin(t);


%Heart
% x = 0.08*(16*sin(t).^3);
% y = 0.08*(13*cos(t) - 5*cos(2*t) - 2*cos(3*t) - cos(4*t));

%Diamond
% a = 0.7;
% x = a*cos(t).^3;
% y = a*sin(t).^3;

%Rose Curve
% a = 0.6;
% petals = 4;
% x= a*cos(petals.*t).*cos(t);
% y= a*cos(petals.*t).*sin(t);

%Hypotrochoid(sometimes playing around messes up, see https://www.desmos.com/calculator/r8bsgvsszo)
% fr = 1.6;
% rr = 0.2;
% d = 0.5;
% x = 0.5*((fr-rr)*cos(t)+d*cos(((fr-rr)/rr)*t));
% y = 0.5*((fr-rr)*sin(t) - d*sin(((fr-rr)/rr)*t));

p = [x', y'];
q = curvspace(p,n);
xx = q(:, 1);
yy = q(:, 2);
c = xx + i*yy;
% c = 0.75*exp(2i*pi*(1:n)/n);
rr = r*ones(size(c));
N = max(0, round(4+.5*log10(r)));

npts = 3*N+2;
circ = exp((1:npts)'*2i*pi/npts);
z = []; for j=1:n
    z=[z;c(j)+rr(j)*circ]; end
A = [0; -ones(size(z))];
zs = [1.5, 1.5i];

%% allows for multiple point charges
rhs_log = 0;
for i = 1:width(zs)
    rhs_log = rhs_log - log(abs(z-zs(i)));
end
rhs = [0; rhs_log];

for j=1:n
    A = [A [1; log(abs(z-c(j)))]];
    for k = 1:N
        zck = (z-c(j)).^(-k);
        A = [A [0; real(zck)] [0;imag(zck)]];
    end
end
disp(size(A));
disp(size(rhs));
X = A\rhs;
disp(size(X));
e=X(1); X(1) = []; %constant voltage on wires
d = X(1:2*N+1:end); X(1:2*N+1:end) = [];
a = X(1:2:end); b=X(2:2:end);


x = linspace(-1.4, 2.2, 120); y = linspace(-1.8, 1.8, 120);
[xx, yy] = meshgrid(x, y); zz=xx+1i*yy; 
%% allows for multiple point charges
uu = 0;
for i=1:width(zs)
    uu=uu+log(abs(zz-zs(i)));
end


for j=1:n
    uu = uu+d(j)*log(abs(zz-c(j)));
    for k=1:N, zck = (zz-c(j)).^(-k); kk=k+(j-1)*N;
        uu = uu+a(kk)*real(zck)+b(kk)*imag(zck); end
end
theta = linspace(0, 2*pi, 1000);
a=0.15; %distance from charge

% x_circ = a*cos(theta) + real(zs);
% y_circ = a*sin(theta) + imag(zs);
% [in,~] = inpolygon(xx,yy,x_circ,y_circ);      

[grad_xx, grad_yy] = gradient(real(uu), 3.6/120, 3.6/120);

% grad_xx(in)=0;grad_yy(in)=0; %filter for too close to point charge
magFX_grid = sqrt(grad_xx.^2 + grad_yy.^2);
% magFX_grid(in)=max(magFX_grid, [], "all"); %filter for too close to point charge

exp_quantiles = quantile(magFX_grid, [0.025, 0.75], "all");
disp(exp_quantiles);
toobig = magFX_grid>exp_quantiles(2);
magFX_grid(toobig) = exp_quantiles(2);

for j=1:n, uu(abs(zz-c(j))<rr(j)) = NaN; end
z = exp(pi*1i*(-50:50)'/50);
for j=1:n, disk = c(j)+rr(j)*z; fill(real(disk), imag(disk), [1 .7 .7])
    hold on, plot(disk, '-r'), end
contour(xx, yy, real(uu), -2:.1:1.2), colormap([0 0 0]), axis([-1.4 2.2 -1.8 1.8])
axis square, plot(real(zs), imag(zs), '.r')
figure;
    %quiver(grad_xx, grad_yy);
    imagesc(magFX_grid);
colorbar;
axis equal