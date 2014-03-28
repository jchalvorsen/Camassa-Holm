close all;
clear all;
tic
x0 = @(x)  2*exp(-abs(x));
%x0 = @(x) sin(pi*x*30);
%x0 = @(x)  2*exp(-abs(x)) + exp(-abs(x - 5));
ref = @(x,t) 2*exp(-abs(x-2*t));
%ref = @(x,t) 2*exp(-abs(x-2*t)) + exp(-abs(x-t -5));


% initialising constants.
xmin = -10;
xmax = 20;
T = 6;
I = 7:1:11;

% All these calculated from above
M = length(I);
N = 2.^I;
H = (xmax-xmin)./N;
K = H/4;
L = floor(T./K);
iterations = sum(L);

wait = waitbar(0,'Launching nuclear missile');
solution = figure;
errornorms = figure;
hold on
color = jet(M);
convergence = zeros(M,1);
progress = 0;
for i = 1:M
    % for readability
    n = N(i);
    h = H(i);
    k = K(i);
    l = L(i);
    
    x = xmin:h:xmax-h;
    u = x0(x)';
    
    % A - operator to go from u to m
    e = ones(n,1);
    A = spdiags([-e/h^2, (1 + 2/h^2) * e, -e/h^2], -1:1, n, n);
    A(n,1) = -1/h^2;
    A(1,n) = -1/h^2;
    
    % B - backwards difference
    B = spdiags([-e/h, e/h], -1:0, n, n);
    B(1,n) = -1/h;

    % C - central difference
    C = spdiags([-e/(2*h), e/(2*h)], [-1, 1], n, n);
    C(1,n) = -1/(2*h);
    C(n,1) = 1/(2*h);
    
    errors = zeros(1,l);
    for j = 1:l
        % updating waitbar
        if i > 1 && floor((L(i-1) + j)*100/iterations) > progress + 5;
            progress = progress + 5;
            waitbar(progress/100)
        end
        % TEST OF NEW FORMAT
        m = A*u;
        mt = m + k*(- B*(m.*u) - m.*(C*u) );
        u = A\mt;
        % Everythin in one operation below
        %u = u + k * (A\((-B*((A*u).*u)) - (A*u).*(C*u)));
        refs = ref(x,j*k)';
        errors(j) = norm(refs-u,2)/(norm(refs,2));
        %errors(j) = 2 - max(u);
    end
    
    convergence(i) = errors(end) - errors(end-1);
    figure(solution)
    hold on
    plot(x,u,'color',color(i,:))
    figure(errornorms)
    hold on
    dots = linspace(0,100,l);
    plot(dots,errors,'color',color(i,:))
end
figure(solution)
legend(cellstr(num2str((N)')),'Location','NorthWest')
plot(x,refs,'k')
delete(wait)
figure(errornorms)
ylabel('relative fault')
xlabel('timesteps')
legend(cellstr(num2str((N)')),'Location','NorthWest')

figure
plot(log(N),log(convergence), 'b*-')
hold on
plot(log(N),log(fliplr(N)) - max(log(N)) + max(log(convergence)),'r*-')
%plot(log(N),log(fliplr(N.^2) - max(log(N.^2)) + max(log(convergence))),'g*-')
legend('Rate of convergence','dy/dx = 1');


toc