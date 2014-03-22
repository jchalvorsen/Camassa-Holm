close all;
clear all;
tic
x0 = @(x)  2*exp(-abs(x));
%x0 = @(x)  2*exp(-abs(2*x + 10));% + -exp(-abs(x - 10));
ref = @(x,t) 2*exp(-abs(x-2*t));



xmin = -10;
xmax = 20;

L = 10000;
T = 6;
dT = T/L;

I = [8, 10, 12, 14];

errors = zeros(length(I),L);
wait = waitbar(0,'Progress');
figure
hold on
col = repmat(['b', 'g', 'r'], 1,ceil(length(I)/3));
for is = 1:length(I)
    n = 2^I(is);
    h = (xmax-xmin)/(n+1);
    x = xmin:h:xmax;
    u = x0(x);
    % testing matrix buildup:
    % A is operator to go from u to m
    e = ones(n+2,1);
    A = spdiags([-e/h^2, (1 + 2/h^2) * e, -e/h^2], -1:1, n+2, n+2);
    A(n+2,1) = -1/h^2;
    A(1,n+2) = -1/h^2;
  
    % B is backwards difference
    B = spdiags([-e/h, e/h], -1:0, n+2, n+2);
    B(1,n+2) = -1/h;
    
    % C is central difference
    C = spdiags([-e/(2*h), e/(2*h)], [-1, 1], n+2, n+2);
    C(1,n+2) = -1/(2*h);
    C(n+2,1) = 1/(2*h);
    
    for i = 2:L
        waitbar(((is-1)*L + i)/(L*length(I)));
        % TEST OF NEW FORMAT
        u = u + dT * (A\((-B*((A*u')'.*u)')' - (A*u')'*C*u')')';
        refs = ref(x,(i-1)*dT);
        errors(is,i) = norm(refs-u,2)/(norm(refs,2));   
    end
    plot(x,u,col(is))
end
plot(x,refs,'k')
delete(wait);
figure
plot(errors')
ylabel('relative fault')
xlabel('timesteps')

toc