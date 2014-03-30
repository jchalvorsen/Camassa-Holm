function [ ] = paper_method( )
close all;
clear all;
x0 = @(x)  sin(pi*x*2/41); %exp(-abs(x));

%a = 40; % period
xmin = -21;
xmax = 20;
n = 1000;
h = (xmax-xmin)/(n+1)
x = xmin:h:xmax;
initial = x0(x)
%initial = [ 0, 0 , initial(3:end-2), 0, 0 ];
% Camassa-Holm:
% mt = -(mu)x - mux
% m = u - uxx


% we go to time T through l steps:
T = 50;
l = 100;

dT = T/l;

% The scheme we're trying to solve is
% mt = -backward(mu) - m*central(u)
% m = u - backward(forward(u))

figure
plot(x,initial,'b')

u = initial;
q = backward(forward(u,h),h);
m = u-q;
figure
plot(x,+q)

figure
plot(x,m)
%axis([xmin, xmax, -0.00001, 0.00001])
%mt = -backward(m.*u,h) - m.*central(u,h)
mt = m + dT * (-backward(m.*u,h) - m.*central(u,h) );
figure
plot(x,mt)
%plot(x,m,'r')
%plot(x,mt,'g')

kappa = log( (1 + 2*n^2 + sqrt(1 + 4*n^2))/(2*n^2))
c = 1/(1 + 2*n^2*(1-exp(-kappa)))

n = n+2;
G = zeros(1,n);
for i = 1:n;
    G(i) = c*(exp(-kappa*i)+exp(kappa*(i-n))/(1-exp(-kappa*n)));
end
size(fft(u))
size(fft(G))
U = ifft(fft(mt).*fft(G));
figure
plot(x,U,'r')

% uexp = zeros(1,n);
% for i = 1:n+2;
%     sums = 0;
%     for k = 1:n+1
%         sums = sums + (exp(-kappa*(i-k)) + exp(kappa*(i-k-(n+2))))*mt(k);
%     end
%     uexp(i) = sums;
% end
% uexp = c/(1-exp(-n*kappa))*uexp;
% figure
% plot(x, uexp, 'cyan')
end 


function [ T ] = generateInitial( x, a)
    d = min(abs(x), abs(a-x));
    T = cosh(d-a/2)/sinh(a/2);
end

function [ T ] = backward( V, h )
   K = [ V(end), V ];
   T = 1/h*(K(2:end) - K(1:end-1));
   %T = [ 0, T(2:end) ];    %neumann boundary
end

function [ T ] = forward( V, h )
   K = [ V, V(1) ];
   T = 1/h*(K(2:end) - K(1:end-1));
   %T = [ T(1:end-1), 0 ];
end

function [ T ] = central( V, h )
   K = [ V(end), V, V(1) ];
   T = 1/(2*h)*(K(3:end) - K(1:end-2));
   %T = [ 0, T(2:end-1), 0 ];
end