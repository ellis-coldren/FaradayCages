clear;clc;
U = 5;
a = .001;
X = -.03:.003:.03;
Y = 0:.003:.03;
[x,y] = meshgrid(X,Y);
u = U + ((U*a^2)./((x.^2 + y.^2).^2)).*(y.^2 - x.^2);
v = ((U*a^2)./((x.^2 + y.^2).^2)).*(-2.*x.*y);
% if y < a^2 - x.^2
%     u = 0; v = 0;
% end
theta = linspace(0, pi, 1000);
x_circ = a*cos(theta);
y_circ = a*sin(theta);
y_circ(end) = 0; 
[in,~] = inpolygon(x,y,x_circ,y_circ); 
u(in)=NaN;
v(in)=NaN;

figure;
    quiver(x, y, u, v)
    hold on 
    plot(x_circ,y_circ)
    fill(x_circ,y_circ,"r")
axis equal
xlim([-.03 .03]);
ylim([0 .03]);