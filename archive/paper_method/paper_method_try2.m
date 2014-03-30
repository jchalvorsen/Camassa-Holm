function [ ] = try2( )
close all;
clear all;
x0 = @(x)  sin(pi*x);%2*exp(-abs(x));% + 2*exp(-abs(2*x + 10));% + -exp(-abs(x - 10));
ref = @(x,t) 2*exp(-abs(x-2*t));
%sin(pi*2*x/20);
tic

xmin = -10;
xmax = 20;

I = 8:2:12;
L = 5000;
errors = zeros(length(I),L);

for is = 1:length(I)
    for ls = 1:length(L)
        n = 2^I(is);
        h = (xmax-xmin)/(n+1);
        x = xmin:h:xmax;
        %initial = generateInitial(x,a) + generateInitial(x,a/2)*2;
        initial = x0(x);
        %initial = [ 0, 0 , initial(3:end-2), 0, 0 ];
        % Camassa-Holm:
        % mt = -(mu)x - mux
        % m = u - uxx
        
        % testing matrix buildup:
        % A is operator to go from u to m
        e = ones(n+2,1);
        A = spdiags([-e/h^2, (1 + 2/h^2) * e, -e/h^2], -1:1, n+2, n+2);
        A(n+2,1) = -1/h^2;
        A(1,n+2) = -1/h^2;
        %full(A)
        
        % B is backwards difference
        B = spdiags([-e/h, e/h], -1:0, n+2, n+2);
        B(1,n+2) = -1/h;
        %full(B)
        
        % C is central difference
        C = spdiags([-e/(2*h), e/(2*h)], [-1, 1], n+2, n+2);
        C(1,n+2) = -1/(2*h);
        C(n+2,1) = 1/(2*h);
        %full(C)
        
        % we go to time T through l steps:
        T = 6;
        l = L(ls);
        
        dT = T/l;
        
        % The scheme we're trying to solve is
        % mt = -backward(mu) - m*central(u);
        % m = u - backward(forward(u))
        
        
        
        %figure
        %plot(x,initial,'b')
        u = initial;
        % kappa = log( (1 + 2*n^2 - sqrt(1 + 4*n^2))/(2*n^2))
        % c = 1/(1 + 2*n^2*(1-exp(-kappa)))
        
        
        %data = zeros(l,n+2);
        %data(1,:) = u;
        %u2 = u;
        %data2 = zeros(l,n+2);
        %data2(1,:) = u2;
        z.wait = waitbar(0,'gone');
        %refs = zeros(l,n+2);
        %refs(1,:) = u;
        error = zeros(l,1);
        figure
        for i = 2:l
            waitbar(i/l);
            %m = (A*u')';
            %m2 = (A*u2')';
            
            %m = u-backward(forward(u,h),h);
            %m2 = u-backward(forward(u,h),h);
            
            %mt = m + dT * (-backward(m.*max(u,0),h) - forward(m.*min(u,0),h) - m.*central(u,h) );
           
            %mt =  m + dT * ((-B*(m.*u)')' - m*C*u');
            %mt2 =  m2 + dT * (-backward(m2.*u2,h) - m2.*central(u2,h));
            %plot(x,u,'r')
            %axis([xmin xmax -0.5 1])
            %u = (A\mt')';
            %u2 = (A\mt2')';
            %hold off
            
            % TEST OF NEW FORMAT
%             m = (A*u')';
%             size(B)
%             size(m)
%             size(m.*u)
%             size(B*(m.*u)' )
%             size((m*(C*u)'))
%             mt = m + dT * ( -B*(m.*u)' - (m*C*u')')
%             mt = mt'
%             u = A\mt;
            %u = u + dT * (A\((-B*(m.*u)')' - m*C*u')')';
            
            
            u = u + dT * (A\((-B*((A*u')'.*u)')' - (A*u')'*C*u')')';
            %q = (A\((-B*((A*u')'.*u)')' - (A*u')'*C*u')')';
            %hold off
            %plot(x,u,'b');
            %hold on
            %plot(x,q,'r');
            %hold off
            % plot reference solution
            refs = ref(x,(i-1)*dT);
            %plot(x,refs,'g')
            %pause(0.001)
            %data(i,:) = u;
            %data2(i,:) = u2;
            errors(is,i) = norm(u-refs,2);
            
        end
        
        delete(z.wait);
        %figure
        %surf(data(1:20:end,:))
        %shading flat
        
        %figure
        %surf(data2(1:20:end,:));
        %shading flat
        
        %if 1-max(u) > 0
        %    errors(is,ls) = 1-max(u)
        %end
        %areaUnderCurveInitial = sum(initial) / length(initial)
        %areaUnderFinalCurve = sum(data(l,:)) / (n+2)
    end
    
    
end

figure
plot(errors')
ylabel('x')
xlabel('timesteps')
shading flat
toc
end 

function [ T ] = generateInitial( x, a)
    d = min(abs(x), abs(a-x));
    T = cosh(d-a/2)/sinh(a/2);
end


function [ T ] = backward( V, h )
   K = [ V(end), V ];
   T = 1/h*(K(2:end) - K(1:end-1));
   %T = [ 0, T(2:end-1), 0 ];    %neumann boundary
end

function [ T ] = forward( V, h )
   K = [ V, V(1) ];
   T = 1/h*(K(2:end) - K(1:end-1));
   %T = [ 0, T(2:end-1), 0 ];    %neumann boundary
end

function [ T ] = central( V, h )
   K = [ V(end), V, V(1) ];
   T = 1/(2*h)*(K(3:end) - K(1:end-2));
   %T = [ 0, T(2:end-1), 0 ];    %neumann boundary
end