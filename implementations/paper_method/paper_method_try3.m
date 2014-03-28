close all;
clear all;
tic
x0 = @(x)  2*exp(-abs(x));
%x0 = @(x) sin(pi*x*30);
%x0 = @(x)  2*exp(-abs(x)) + exp(-abs(x - 5));
ref = @(x,t) 2*exp(-abs(x-2*t));
%ref = @(x,t) 2*exp(-abs(x-2*t)) + exp(-abs(x-t -5));

%m0 = @(x) sech(x).^2;

xmin = -10;
xmax = 20;

L = 1000;
T = 4;
dT = T/L;

I = 6:1:14;

%errors = zeros(L,length(I));
wait = waitbar(0,'Progress');
sol = figure;
fault = figure;
hold on
col = jet(length(I));
conv = zeros(length(I),1);
for is = 1:length(I)
    n = 2^I(is);
    h = (xmax-xmin)/(n+1);
    dT = h/8;
    L = floor(T/dT);
    x = xmin:h:xmax;
    u = x0(x);
    % testing matrix buildup:
    % A is operator to go from u to m
    e = ones(n+2,1);
    A = spdiags([-e/h^2, (1 + 2/h^2) * e, -e/h^2], -1:1, n+2, n+2);
    A(n+2,1) = -1/h^2;
    A(1,n+2) = -1/h^2;
    %u = (A\m0(x)')';
    % B is backwards difference
    B = spdiags([-e/h, e/h], -1:0, n+2, n+2);
    B(1,n+2) = -1/h;

    % C is central difference
    C = spdiags([-e/(2*h), e/(2*h)], [-1, 1], n+2, n+2);
    C(1,n+2) = -1/(2*h);
    C(n+2,1) = 1/(2*h);
    
    errors = zeros(1,L);
    for i = 2:L
        if mod(((is-1)*L + i)/(L*length(I))*100,5) == 0
            waitbar(((is-1)*L + i)/(L*length(I)));
        end
        % TEST OF NEW FORMAT 
        u = u + dT * (A\((-B*((A*u')'.*u)')' - (A*u')'*C*u')')';
        refs = ref(x,(i-1)*dT);
        errors(i) = norm(refs-u,2)/(norm(refs,2));   
    end
    conv(is) = errors(end)-errors(end-1);
    figure(sol)
    hold on
    plot(x,u,'color',col(is,:))
    figure(fault)
    hold on
    dots = linspace(0,100,L);
    plot(dots,errors,'color',col(is,:))
end
figure(sol)
legend(cellstr(num2str((2.^I)')),'Location','NorthWest')
plot(x,refs,'k')
delete(wait)
figure(fault)
ylabel('relative fault')
xlabel('timesteps')
legend(cellstr(num2str((2.^I)')),'Location','NorthWest')

figure
loglog(2.^I,conv, 'b*-');
hold on
loglog(2.^I,fliplr(2.^(I))/100,'r*-')

toc