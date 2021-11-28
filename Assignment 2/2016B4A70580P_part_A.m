clear all
format long
clc

N = 200; % no of discritization points
mm = 4;
[x1,Y1] = solve(0.068,200);
[x2,Y2] = solve(0.1,200);
[x3,Y3] = solve(0.5,200);

y1 = zeros(N+1,1);
y2 = zeros(N+1,1);
y3 = zeros(N+1,1);
y4 = zeros(N+1,1);
for i = 1:N+1
    y1(i) = Y1(mm*(i-1)+1);
    y2(i) = Y1(mm*(i-1)+2);
    y3(i) = Y1(mm*(i-1)+3);
    y4(i) = Y1(mm*(i-1)+4);
end
hold on

xlabel('$y$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex');
ylabel('$q_y(y)$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex');
title('$q_y(y)$ Plot for Different values of $Kn$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex');

col1 = plot(x1,y4, 'Color','b');

for i = 1:N+1
    y1(i) = Y2(mm*(i-1)+1);
    y2(i) = Y2(mm*(i-1)+2);
    y3(i) = Y2(mm*(i-1)+3);
    y4(i) = Y2(mm*(i-1)+4);
end
col2 = plot(x2,y4, 'Color','g');

for i = 1:N+1
    y1(i) = Y3(mm*(i-1)+1);
    y2(i) = Y3(mm*(i-1)+2);
    y3(i) = Y3(mm*(i-1)+3);
    y4(i) = Y3(mm*(i-1)+4);
end
col3 = plot(x3,y4, 'Color','r');
legend([col1(1),col2(1),col3(1)], 'Kn=0.068', 'Kn=0.1', 'Kn=0.5', 'Location','northeast');

function [x,y] = solve(Kval,N)

    Kn = Kval;
    F = 0.23;
    eps = 0.000000005;

%     N = 200; % no of discritization points
    mm = 4;
    h = 1/N; % mesh size
    x = (0:h:1)'; % nodes
    M = zeros(mm*(N+1),mm*(N+1));
    b = zeros(mm*(N+1),1); %M*y=b;

    Done = 0;
    Max_iter = 10;

    y0 = zeros(mm*(N+1),1); %Initial guess
    for i = 1:N+1
       y0(mm*(i-1)+1) = 0;
       y0(mm*(i-1)+2) = 0;
       y0(mm*(i-1)+3) = 1;
       y0(mm*(i-1)+4) = 0;
    end

    iter = 0;

    while (~Done)
        for i = 1:N
            xip1 = 0.5*(x(i)+x(i+1));
            yip1 = 0.5*(y0(mm*(i-1)+1:mm*i)+y0(mm*i+1:mm*(i+1)));
            Aip1 = [0, -1/Kn, 0, 0; 0, 0, -F/(yip1(3)^2+eps), 0; 0, 0, 0, -4/(15*Kn); 0, 2*yip1(2)/Kn, 0, 0];
            Qip1 = [-yip1(2)/Kn; F/(yip1(3)+eps); (-4/(15*Kn))*yip1(4); yip1(2)*yip1(2)/Kn];
            fip1 = Qip1-Aip1*yip1;
            M(mm*(i-1)+1:mm*i,mm*(i-1)+1:mm*i) = -(0.5*Aip1+eye(mm)/h);
            M(mm*(i-1)+1:mm*i,mm*i+1:mm*(i+1)) = (-0.5*Aip1+eye(mm)/h);
            b(mm*(i-1)+1:mm*i) = fip1;
        end 

        ya = y0(1:mm);
        yb = y0(mm*N+1:mm*(N+1));
        Ba = [(2/(pi*(ya(3)+eps)))^0.5, 1, -0.5*ya(1)*(2/(pi*(ya(3)^3+eps)))^0.5, 0;
        0, 0, ((2/pi)^0.5)*((1/(ya(3)+eps))^0.5+(1/(ya(3)^3+eps))^0.5), 1;
        0, 0, 0, 0;
        0, 0, 0, 0];
        Bb = [0, 0, 0, 0;
        0, 0, 0, 0;
        -1*(2/(pi*(yb(3)+eps)))^0.5, 1, 0.5*yb(1)*(2/(pi*(yb(3)^3+eps)))^0.5, 0;
        0, 0, -1*((2/pi)^0.5)*((1/(ya(3)+eps))^0.5+(1/(ya(3)^3+eps))^0.5), 1];

        M(mm*((N+1)-1)+1:mm*(N+1),mm*((N+1)-1)+1:mm*(N+1)) = Bb;
        M(mm*((N+1)-1)+1:mm*(N+1),1:mm) = Ba;
        alpha = Ba*ya+Bb*yb-[ya(2)+ya(1)*(2/(pi*(ya(3)+eps)))^0.5; ya(4)+2*(ya(3)-1)*(2/(pi*(ya(3)+eps)))^0.5; yb(2)-yb(1)*(2/(pi*(yb(3)+eps)))^0.5; yb(4)-2*(yb(3)-1)*(2/(pi*(yb(3)+eps)))^0.5];
        b(mm*((N+1)-1)+1:mm*(N+1)) = alpha;

        y = M\b;

        error = max(abs(y-y0));% max |y_i-y_i^*|
        iter = iter+1;
        disp(['Error in step ',num2str(iter),' is : ',num2str(error)])
        y0 = y;

        Done = (iter>=Max_iter)||(error<0.000001);
    end
end
