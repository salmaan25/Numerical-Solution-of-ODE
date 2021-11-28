% ASR NODES 2020
% The non-linear mid-point method
% SOLVES y'(x)=q(x,y) 
% g(y(a),g(b)) = 0 

clear all
format long
clc

N = 50; % no of discritization points
mm = 2;
h = 1/N; % mesh size
u0=1; % u at left boundary
u1=1; % u at right boundary
x = (0:h:1)'; % nodes
M = zeros(mm*(N+1),mm*(N+1));
b = zeros(mm*(N+1),1); %M*y=b;

Done = 0;
Max_iter = 10;

omega =-0.5;
y0 = zeros(mm*(N+1),1); %Initial guess
for i = 1:N+1
   y0(mm*(i-1)+1) = 1;
   y0(mm*(i-1)+2) = 0;
end

iter = 0;

while (~Done)
    for i = 1:N
        xip1 = 0.5*(x(i)+x(i+1));
        yip1 = 0.5*(y0(mm*(i-1)+1:mm*i)+y0(mm*i+1:mm*(i+1)));
        Tip1 = yip1(1);
        qip1 = yip1(2);
        Aip1 = [-omega*qip1*Tip1^(omega-1), -Tip1^(omega-1);
            0, 0];
        Qip1 = [-qip1*Tip1^omega; sin(2*pi*xip1)];
        fip1 = Qip1-Aip1*yip1;
        M(mm*(i-1)+1:mm*i,mm*(i-1)+1:mm*i) = -(0.5*Aip1+eye(mm)/h);
        M(mm*(i-1)+1:mm*i,mm*i+1:mm*(i+1)) = (-0.5*Aip1+eye(mm)/h);
        b(mm*(i-1)+1:mm*i) = fip1;
    end 
    
    Ba = [1, 0;
    0, 0];
    Bb = [0, 0;
    1, 0];
     
    M(mm*((N+1)-1)+1:mm*(N+1),mm*((N+1)-1)+1:mm*(N+1)) = Bb;
    M(mm*((N+1)-1)+1:mm*(N+1),1:mm) = Ba;
    alpha = [u0; u1];
    b(mm*((N+1)-1)+1:mm*(N+1)) = alpha;
    
    y = M\b;
    
    error = max(abs(y-y0));% max |y_i-y_i^*|
    iter = iter+1;
    disp(['Error in step ',num2str(iter),' is : ',num2str(error)])
    y0 = y;

    Done = (iter>=Max_iter)||(error<0.000001);
    
    
end


figure(1);
y1 = zeros(N+1,1);
y2 = zeros(N+1,1);
for i = 1:N+1
    y1(i) = y(mm*(i-1)+1);
    y2(i) = y(mm*(i-1)+2);
end

plot(x,y1,'ok','LineWidth',1),grid on;
xlabel('$x$','FontSize',13,'Color','k', 'Interpreter', 'latex')
ylabel('$T(x)$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex')
