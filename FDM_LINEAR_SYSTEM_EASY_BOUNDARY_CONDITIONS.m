% ASR NODES 2020
% SOLVES y'(t)=A(t)y+f(t) 


format long
N = 50; % no of discritization points
mm = 2;
h = 1/N; % mesh size
u0=0; % u at left boundary
u1=0; % u at right boundary
x = (h:h:1-h)'; % interior nodes
M = zeros(mm*(N-1),mm*(N-1));
b = zeros(mm*(N-1),1);

%  left
A = [0, 1;
    0, 0];
B0 = [0, 0;
    0, 1];
al0 = [u0; 0];
M(1:mm,1:mm) = -(A+B0/h);
M(1:mm,mm+1:mm*2) = (eye(mm)+B0)/(2*h);
b(1:mm) = al0/(2*h)+[0; exp(x(1))];
for i = 2:N-2
    A = [0, 1;
        0, 0];
    M(mm*(i-1)+1:mm*i,mm*(i-1)+1:mm*i) = -A;
    M(mm*(i-1)+1:mm*i,mm*(i-2)+1:mm*(i-1)) = -eye(mm)/(2*h);
    M(mm*(i-1)+1:mm*i,mm*i+1:mm*(i+1)) = eye(mm)/(2*h);
    b(mm*(i-1)+1:mm*i) = [0; exp(x(i))];
end
A = [0, 1;
        0, 0];
B1 = [1, 0;
    0, 0];
al1 = [0; u1];
M(mm*((N-1)-1)+1:mm*(N-1),mm*((N-1)-1)+1:mm*(N-1)) = -(A-B1/h);
M(mm*((N-1)-1)+1:mm*(N-1),mm*((N-1)-2)+1:mm*((N-1)-1)) = -(eye(mm)+B1)/(2*h);
b(mm*((N-1)-1)+1:mm*(N-1)) = [0; exp(x(N-1))]-al1/(2*h);

u = M\b;


figure(1);
y1 = zeros(N-1,1);
y2 = zeros(N-1,1);
for i = 1:N-1
    y1(i) = u(mm*(i-1)+1);
    y2(i) = u(mm*(i-1)+2);
end

plot(x,y1,'ok','LineWidth',1),grid on;
hold on
plot(x,-1+exp(x)+u0-exp(1)*x+u1*x,'r','LineWidth',1.5)