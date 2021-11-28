format long
mesh = [12,25,50,100,200,400];
Error = zeros(size(mesh));

for i = 1:length(mesh)
    Error(i) =  FDM261020OoC(mesh(i));
end

figure(1);
loglog(mesh,Error,'-ok','LineWidth',1)
grid on;

xlabel('$N$','FontSize',13,'Color','k', 'Interpreter', 'latex')
ylabel('$\epsilon_{h}$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex')
title('error in log-log plot','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex')


figure(2);
plot(mesh,Error.*((mesh.^2)),'-ok','LineWidth',1)
grid on;



function error = FDM261020OoC(N)
h = 1/N;
u0=0;
u1=0;
x = (h:h:1-h)';
M = zeros(N-1,N-1);
b = zeros(N-1,1);

M(1, 1) = -2/(h^2);
M(1, 2) = 1/(h^2);
b(1) = exp(x(1))-u0/(h^2);

for i = 2:N-2
    M(i, i-1) = 1/(h^2);
    M(i, i+1) = 1/(h^2);
    M(i, i) = -2/(h^2);
    b(i) = exp(x(i));
end

M(N-1, N-1) = -2/(h^2);
M(N-1, N-2) = 1/(h^2);
b(N-1) = exp(x(N-1))-u1/(h^2);

u = M\b;

error = max(abs(u-(-1+exp(x)+u0+x-exp(1)*x-u0*x+u1*x)));

end