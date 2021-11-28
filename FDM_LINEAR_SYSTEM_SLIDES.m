% ASR NODES 2020
% The linear mid-point method
% SOLVES y'(x)=A(x)y+f(x) 
% Ba*y0+Bb*yN = alpha 


format long
N = 50; % no of discritization points
mm = 2;
h = 1/N; % mesh size
u0=1; % u at left boundary
u1=1; % u at right boundary
x = (0:h:1)'; % nodes
M = zeros(mm*(N+1),mm*(N+1));
b = zeros(mm*(N+1),1);

%  left
A = [0, -1;
    0, 0];
Ba = [1, 0;
    0, 0];
Bb = [0, 0;
    1, 0];
alpha = [u0, u1]';

for i = 1:N
    xip1 = 0.5*(x(i)+x(i+1));
    M(mm*(i-1)+1:mm*i,mm*(i-1)+1:mm*i) = -(0.5*A+eye(mm)/h);
    M(mm*(i-1)+1:mm*i,mm*i+1:mm*(i+1)) = (-0.5*A+eye(mm)/h);
    b(mm*(i-1)+1:mm*i) = [0; sin(2*pi*xip1)];
end

M(mm*((N+1)-1)+1:mm*(N+1),mm*((N+1)-1)+1:mm*(N+1)) = Bb;
M(mm*((N+1)-1)+1:mm*(N+1),1:mm) = Ba;
b(mm*((N+1)-1)+1:mm*(N+1)) = alpha;


u = M\b;


figure(1);
y1 = zeros(N+1,1);
y2 = zeros(N+1,1);
for i = 1:N+1
    y1(i) = u(mm*(i-1)+1);
    y2(i) = u(mm*(i-1)+2);
end

plot(x,y1,'ok','LineWidth',1),grid on;
hold on

plot(x,(4*pi^2+sin(2*pi*x))/(4*pi^2),'r','LineWidth',1.5)
xlabel('$x$','FontSize',13,'Color','k', 'Interpreter', 'latex')
ylabel('$T(x)$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex')
max(abs(y1-(4*pi^2+sin(2*pi*x))/(4*pi^2)))