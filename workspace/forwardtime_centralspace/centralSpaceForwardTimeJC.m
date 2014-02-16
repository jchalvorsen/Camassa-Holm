function [ ] = centralSpaceForwardTime( )
x0 = @(x) 2*exp(-abs(x)) + 5*exp(-abs(x-60));
%x0 = @(x) 1*exp(-(x-2).^2) + 4*exp(-(x-50).^2);
x0 = @(x) 10./(3+abs(x)).^2;

close all
n = 2^10;
xmax = 160;
xmin = -50;
%x0 = @(x) sin(x*pi/(xmax-xmin));
h = (xmax-xmin)/(n-1);
x = xmin:h:xmax;
initial = x0(x);
plot(x,initial)



k = 0.1;
A = buildMatrix(n, h);
%full(A)

T = 1500;      % number of iterations for solutions
% finding F
F = zeros(T,n);
U = zeros(T,n);
prev = initial;
for i = 1:T
    c1 = firstCentral(prev, h);
    c2 = secondCentral(prev, h);
    c3 = thirdCentral(prev, h);
    F(i,:) = c2.*(k + 2*k*c1) + prev.*(1+k*c3-3*k*c1);
    % Finding U
    U(i,:) = A\F(i,:)';
    U(i,:) = U(i,:) - min(U(i,:));
    U(i,:) = U(i,:)./(max(U(i,:)));
    prev = U(i,:);
end
U = [initial; U];
k = figure;
axis manual
axis([ xmin xmax 0 1])

%hold on
for i = 1:T
    plot(x,U(i,:))
    axis([ xmin xmax 0 1])
    pause(0.01)    
end
%figure
%surf(U)
%shading flat
%axis manual
%surf(U)
end

function [ T ] = firstCentral( V, h )

K = [ V(end), V, V(1) ];
T = 1/(2*h)*(K(3:end) - K(1:end-2));
%T = (K(3:end) - K(1:end-2));

end

function [ T ] = secondCentral( V, h )

K = [ V(end), V, V(1) ];
T = 1/h^2*( K(3:end) - 2*K(2:end-1) + K(1:end-2));
%T = ( K(3:end) - 2*K(2:end-1) + K(1:end-2));

end

function [ T ] = thirdCentral( V, h )

K = [ V(end-1:end), V, V(1:2) ];
T = 1/(2*h^2)* ( K(5:end) - 2*K(4:end-1) + 2*K(2:end-3) - K(1:end-4));

end

function [ A ] = buildMatrix( n, h )
e = ones(n,1);
A = spdiags([-1/h^2*e, (1+2/h^2)*e, -1/h^2*e], -1:1, n, n);
A(1,n) = -1/h^2;
A(n,1) = -1/h^2;
end