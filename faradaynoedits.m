n = 12;r = 0.1;
c = exp(2i*pi*(1:n)/n);
rr = r*ones(size(c));
N = max(0, round(4+.5*log10(r)));

npts = 3*N+2;
circ = exp((1:npts)'*2i*pi/npts);
z = []; for j=1:n
    z=[z;c(j)+rr(j)*circ]; end
disp(size(c));
disp(size(z));
A = [0; -ones(size(z))];
zs = 1.5;
rhs = [0; -log(abs(z-zs))];
disp(size(log(abs(z-c(1)))));
disp(size(rhs));
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
[xx, yy] = meshgrid(x, y); zz=xx+1i*yy; uu=log(abs(zz-zs));
% disp(uu)

for j=1:n
    uu = uu+d(j)*log(abs(zz-c(j)));
    for k=1:N, zck = (zz-c(j)).^(-k); kk=k+(j-1)*N;
        uu = uu+a(kk)*real(zck)+b(kk)*imag(zck); end
end

for j=1:n, uu(abs(zz-c(j))<rr(j)) = NaN; end
z = exp(pi*1i*(-50:50)'/50);
for j=1:n, disk = c(j)+rr(j)*z; fill(real(disk), imag(disk), [1 .7 .7])
    hold on, plot(disk, '-r'), end
contour(xx, yy, real(uu), -2:.1:1.2), colormap([0 0 0]), axis([-1.4 2.2 -1.8 1.8])
axis square, plot(real(zs), imag(zs), '.r')
