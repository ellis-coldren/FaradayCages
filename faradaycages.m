n = 30; %number of wires
r = 0.02; %radius of wires
t = linspace(0, 2*pi, 1000); %setting parameter t from 0 to 2pi
using_vec = 1;

%Circle
% radius = 0.75;
% x = radius*cos(t);
% y = radius*sin(t);

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

% -------------------- Reading .vec files ----------------------- %
vecFile = 'fruit.vec';
vecRead = xmlread(vecFile);
using_vec = 0;

sketch_edges = vecRead.getElementsByTagName('edge');
x_shape = [];
y_shape = [];

for i = 0:sketch_edges.getLength-1
    edge = sketch_edges.item(i);
    % the xywdense has format: xywdense(layer(int) x1, y1, w1, x2, y2, w2, ...);
    densesample = char(edge.getAttribute('curve'));
    sample_substring = regexp(densesample, '[+-]?\d*\.?\d+', 'match');
    sample_double = str2double(sample_substring);

    length_x = length(sample_double(2:3:end));
    length_y = length(sample_double(3:3:end));

    x_indices = round(linspace(1, length_x, n));
    y_indices = round(linspace(1, length_y, n));
    
    x_sample = sample_double(2:3:end);
    x_to_add = x_sample(x_indices);
    y_sample = sample_double(3:3:end);
    y_to_add = y_sample(y_indices);

    x_shape = [x_shape, x_to_add];
    y_shape = [y_shape, y_to_add];
end
min_x = min(x_shape);
max_x = max(x_shape);
min_y = min(y_shape);
max_y = max(y_shape);
x_scale = (3 - (-3))/(max_x - min_x);
y_scale = (3 - (-3))/(max_y - min_y);
x = x_scale*x_shape;
y = -y_scale*y_shape;

% --------------------------------------------------------------- %

% p = [x', y'];
% q = curvspace(p,n+1);
% xx = q(:, 1);
% yy = q(:, 2);
% c = xx + 1i*yy;

p = [x', y'];
q = curvspace(p,n+1); % generates points that interpolate curve
q = unique(q, 'rows')
xx = q(:, 1);
yy = q(:, 2);
if using_vec == 0
    disp('using vec');
    xx = p(:, 1);
    yy = p(:, 2);
end
c = xx + 1i*yy; % curve as vector of complex coordinates

% c = 0.75*exp(2i*pi*(1:n)/n);
rr = r*ones(size(c));
N = max(0, round(4+.5*log10(r)));

npts = 3*N+2;
circ = exp((1:npts)'*2i*pi/npts);
z = []; for j=1:n
    z=[z;c(j)+rr(j)*circ]; end
A = [0; -ones(size(z))];

figure;
hold on;
% plot(real(c), imag(c), 'k-', 'LineWidth', 2); % Plot the main curve
plot(real(z), imag(z), 'b.'); % Plot the wires
axis equal;
axis([-1.8 1.8 -1.8 1.8]);
grid on;
title('Curve and Wires - Click to Select Points');

disp('Click on the plot to add points. Press Enter to stop.');
[x_click, y_click] = ginput; % Capture clicked points

% Convert clicked points to complex numbers
clicked_points = x_click + 1i*y_click;

% Append clicked points to the zs list
zs = [clicked_points'];
disp(zs);
% allows for multiple point charges
rhs_log = 0;
for j = 1:width(zs)
    rhs_log = rhs_log - ((-1).^j)*log(abs(z-zs(j)));
end
rhs = [0; rhs_log];

for j=1:n
    A = [A [1; log(abs(z-c(j)))]];
    for k = 1:N
        zck = (z-c(j)).^(-k);
        A = [A [0; real(zck)] [0;imag(zck)]];
    end
end

X = A\rhs;
e=X(1); X(1) = []; %constant voltage on wires
d = X(1:2*N+1:end); X(1:2*N+1:end) = [];
a = X(1:2:end); b=X(2:2:end);


min_x = min(x);
max_x = max(x);
min_y = min(y);
max_y = max(y);


x = linspace(min_x-2, max_x+2, 120); y = linspace(min_y-2, max_y+2, 120);
[xx, yy] = meshgrid(x, y); zz=xx+1i*yy; 
% allows for multiple point charges
uu = 0;
for j=1:width(zs)
    uu=uu+(-1).^j*log(abs(zz-zs(j)));
end


for j=1:n
    uu = uu+d(j)*log(abs(zz-c(j)));
    for k=1:N, zck = (zz-c(j)).^(-k); kk=k+(j-1)*N;
        uu = uu+a(kk)*real(zck)+b(kk)*imag(zck); end
end

% theta = linspace(0, 2*pi, 1000);
% a=0.15; %distance from charge
% x_circ = a*cos(theta) + real(zs);
% y_circ = a*sin(theta) + imag(zs);
% [in,~] = inpolygon(xx,yy,x_circ,y_circ);      

[grad_xx, grad_yy] = gradient(real(uu), abs(max_x-min_x)/120, abs(max_y-min_y)/120);


%grad_xx(in)=0;grad_yy(in)=0; %filter for too close to point charge
magFX_grid = sqrt(grad_xx.^2 + grad_yy.^2);
%magFX_grid(in)=max(magFX_grid, [], "all"); %filter for too close to point charge

exp_quantiles = quantile(magFX_grid, [0.025, 0.65], "all");
toobig = magFX_grid>exp_quantiles(2);
magFX_grid(toobig) = exp_quantiles(2);

for j=1:n, uu(abs(zz-c(j))<rr(j)) = NaN; end
z = exp(pi*1i*(-50:50)'/50);
for j=1:n, disk = c(j)+rr(j)*z; fill(real(disk), imag(disk), [1 .7 .7])
    hold on, plot(disk, '-r'), end
contour(xx, yy, real(uu), -2:.1:2), colormap([0 0 0]), axis([min_x max_x min_y max_y])
axis square, plot(real(zs), imag(zs), '.r')
%, plot(real(zs(2)), imag(zs(2)), '.r')
figure;
    %quiver(grad_xx, grad_yy);
    imagesc(magFX_grid);
    set(gca, 'YDir', 'normal'); 
colorbar;
axis equal