clear all
format long
clc

mm = 4;
Nf = 800;
[xf, yf] = solve(0.068, Nf);
y1f = zeros(1,Nf+1); %Check this with zeros(N+1,1)
for j = 1:Nf+1
    y1f(j) = yf(mm*(j-1)+1);
end

mesh = [25,50,100,200]';
ratio =2*[16,8,4,2]';
Error1 = zeros(4,1);


for i = 1:4
    [x, y] = solve(0.068,mesh(i));
    y1 = zeros(1,mesh(i)+1); %Check this with zeros(N+1,1)
    for j = 1:mesh(i)+1
        y1(j) = y(mm*(j-1)+1);
    end
    k = ratio(i);
    Error1(i) = max(abs(y1f(1:k:Nf+1)-y1)); % l_inf error in y1
end

% hold on
p1 = polyfit(log(1./mesh),log(Error1(:,1)),1);
% title(['$N$ vs $Error$'],'FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex');
xlabel('$N$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex');
ylabel('$Error$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex');
title(['$\log(N)$ vs $\log(Error)$, $Kn = $', num2str(0.068), ', Slope of best fit line = ', num2str(p1(1))],'FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex');
col1 = loglog(mesh,Error1(:,1),'-ok','Color','g');
legend([col1(1)], 'Kn=0.068', 'Location','northeast');

[xf, yf] = solve(0.5, Nf);
y1f = zeros(1,Nf+1); %Check this with zeros(N+1,1)
for j = 1:Nf+1
    y1f(j) = yf(mm*(j-1)+1);
end
Error2 = zeros(4,1);
for i = 1:4
    [x, y] = solve(0.5,mesh(i));
    y1 = zeros(1,mesh(i)+1); %Check this with zeros(N+1,1)
    for j = 1:mesh(i)+1
        y1(j) = y(mm*(j-1)+1);
    end
    k = ratio(i);
    Error2(i) = max(abs(y1f(1:k:Nf+1)-y1)); % l_inf error in y1
end

p2 = polyfit(log(1./mesh),log(Error2(:,1)),1);
% title(['$\log(N)$ vs $\log(Error)$, $Kn = $', num2str(0.5), ', Slope of best fit line = ', num2str(p2(1))],'FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex');
% col2 = loglog(mesh,Error2(:,1),'-ok','Color','r');
% legend([col2(1)], 'Kn=0.5', 'Location','northeast');

% legend([col1(1),col2(1)], 'Kn=0.068', 'Kn=0.5', 'Location','northeast');


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
