% ASR NODES 2020
% SOLVES y''(t)=exp(x), y(0)=u0, y'(1)=u1, i.e.,
% Dirichlet boundary condition on the left-side of the domain
% Neumann boundary condition  on the right-side of the domain
% extrapolation: u'[x] = (1.5*u[x]-2*u[x-h]+0.5*u[x-2*h])/h;

format long
mesh = [12,25,50,100,200,400];
Error = zeros(size(mesh));

for i = 1:length(mesh)
    Error(i) =  FDM261020OoC(mesh(i));
end


p = polyfit(log(1./mesh),log(Error),1);
 
disp(['Order of the method is: ', num2str(p(1))])

figure(1);
loglog(mesh,Error,'-ok','LineWidth',1)
grid on;

xlabel('$N$','FontSize',13,'Color','k', 'Interpreter', 'latex')
ylabel('$\epsilon_{h}$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex')
title(['Order of the method is: ', num2str(p(1))],'FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex')




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

M(N-1, N-1) = -2/(3*(h^2));
M(N-1, N-2) =  2/(3*(h^2));
b(N-1) = exp(x(N-1))-2*u1/(3*h);

u = M\b;

error = max(abs(u-(-1+exp(x)+u0-exp(1)*x+u1*x)));

end
