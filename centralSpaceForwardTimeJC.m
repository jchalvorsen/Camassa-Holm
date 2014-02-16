function [ ] = centralSpaceForwardTime( )
x0 = @(x) 2*exp(-abs(x-2)) + exp(-abs(x-5));
close all
n = 500;
xmax = 20;
xmin = -10;
x = linspace(xmin, xmax, n)
initial = x0(x);
plot(x,initial)
hold on
T = thirdCentral(initial,1/n^(1/2))
plot(x,T)


end

function [ T ] = firstCentral( V, h )

K = [ V(end), V, V(1) ];
T = 1/(2*h)*(K(3:end) - K(1:end-2));

end

function [ T ] = secondCentral( V, h )

K = [ V(end), V, V(1) ];
T = 1/h^2*( K(3:end) - 2*K(2:end-1) + K(1:end-2));

end

function [ T ] = thirdCentral( V, h )

K = [ V(end-1:end), V, V(1:2) ];
T = 1/(2*h^2)* ( K(5:end) - 2*K(4:end-1) + 2*K(2:end-3) - K(1:end-4));

end