clear all
close all
c = 1;

u = @(x,t) exp(-abs(x-c*t));
u2 = @(x,t) exp(-abs(+10 + x-2*c*t));

n = 200;
x = linspace(0,25,n);
t = linspace(0,25,n);
U = zeros(length(x),length(t));


for i = 1:n
    for j = 1:n
        
        temp = u(x(i),t(j)) + u2(x(i),t(j));
        U(i,j) = temp;
        %plot(U(i,j));
        %pause(0.1)
    end
    
end

for i = 1:n
    axis([0 2 0 n])
    plot(U(:,i))
    pause(0.1)
end

figure
shading flat;
surf(x,t,U)
view(140, 30)
set(gca,'FontSize',15)
xlabel('space x')
ylabel('time t')
zlabel('u')