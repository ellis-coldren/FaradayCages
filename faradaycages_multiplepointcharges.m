n = 40;
r = 0.01;
c = exp(2i*pi*(1:n)/n);
rr = r*ones(size(c));
N = max(0, round(4+.5*log10(r)));

npts = 3*N+2;
circ = exp((1:npts)'*2i*pi/npts);
z = []; for j=1:n
    z=[z;c(j)+rr(j)*circ]; end
A = [0; -ones(size(z))];
zs = [1.5, -1.5i];
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

X = A\rhs;

e=X(1); X(1) = []; %constant voltage on wires
d = X(1:2*N+1:end); X(1:2*N+1:end) = [];
a = X(1:2:end); b=X(2:2:end);


x = linspace(-1.4, 2.2, 120); y = linspace(-1.8, 1.8, 120);
[xx, yy] = meshgrid(x, y); zz=xx+1i*yy; 
uu = 0;
for i=1:width(zs)
    uu=uu+log(abs(zz-zs(i)));
end

for j=1:n
    uu = uu+d(j)*log(abs(zz-c(j)));
    for k=1:N, zck = (zz-c(j)).^(-k); kk=k+(j-1)*N;
        uu = uu+a(kk)*real(zck)+b(kk)*imag(zck); end
end


[grad_xx, grad_yy] = gradient(real(uu), xx(1, :), yy(:, 1));
magFX_grid = sqrt(grad_xx.^2 + grad_yy.^2);

for j=1:n, uu(abs(zz-c(j))<rr(j)) = NaN; end
z = exp(pi*1i*(-50:50)'/50);
for j=1:n, disk = c(j)+rr(j)*z; fill(real(disk), imag(disk), [1 .7 .7])
    hold on, plot(disk, '-r'), end
contour(xx, yy, real(uu), -2:.1:2), colormap([0 0 0]), axis([-1.4 2.2 -1.8 1.8])
axis square, plot(real(zs(1)), imag(zs(1)), '.r'), plot(real(zs(2)), imag(zs(2)), '.r')