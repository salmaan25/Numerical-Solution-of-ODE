function [] = shooting_method_p2()

    clc

    global Re B0 Bb eps;
%     N = 100;
    eps = 0.0001;
    tspan=[0 1]; % set time interval
    Re = 10;%Reynolds number 

    B0 = [1, 0, 0, 0;
        0, 1, 0, 0;
        0, 0, 0, 0;
        0, 0, 0, 0]; % B0 boundary matrix at t = 0;
    Bb = [0, 0, 0, 0;
        0, 0, 0, 0;
        1, 0, 0, 0;
        0, 1, 0, 0];  % Bb boundary matrix at t = b;
    b1 = 0;
    b2 = 0;
    b3=1;
    b4=0;
    alpha = [b1 b2 b3 b4]';

    I = eye(4);

    figure1 = figure('Units','inches','Position',[0,0,6,6*0.7]); 

    % first shooting
    s0  = [b1 b2 b3 b4]'; %inital guss for the y(0) where y is the solution
    hold on
    for i = 1:10
        [t,y1]=classicalRK_4(@combofun,tspan, [s0; I(:,1)], 100);
        [t,y2]=classicalRK_4(@combofun,tspan, [s0; I(:,2)], 100);
        [t,y3]=classicalRK_4(@combofun,tspan, [s0; I(:,3)], 100);
        [t,y4]=classicalRK_4(@combofun,tspan, [s0; I(:,4)], 100);
%         residual = B0*[y1(1,1); y1(2,1); y1(3,1); y1(4,1)]]+Bb*[y1(1,N+1); y1(2,N+1); y1(3,N+1); y1(4,N+1)]-alpha; % The boundary condition
        residual = B0*[y1(1,1); y1(2,1); y1(3,1); y1(4,1)]+Bb*[y3(1,101); y3(2,101); y3(3,101); y3(4,101)]-alpha; % The boundary condition
        zb = [y1(5:8,101), y2(5:8,101), y3(5:8,101), y4(5:8,101)]; % The value of zi at b (the end point)
        snew = s0-(B0+Bb*zb)\residual;
        disp(snew);
%         x=norm(snew-s0);
%         disp("x is ");
%         disp(x);
%         disp("i = ")
%         disp(i);
% 
%         if(abs(x)<=eps)
%             break;
%         end
       s0 = snew;
    end
%     p1 = plot(t,y1(1,:),'g-.');
%     p2=plot(t,y1(2,:),'g-.');
%     p3=plot(t,y1(3,:),'g-.');
    p4=plot(t,y1(4,:),'g-.');
    xlabel('$t$','FontSize',13,'Color','k', 'Interpreter', 'latex')
    ylabel('$y_1(t)$','FontSize',13,'Color','k', 'Interpreter', 'latex')
    title('Shooting method','FontSize',13,'Color','k', 'Interpreter', 'latex')
    % h1 = legend([p1 p2 p0],{'1st shooting','2nd shooting','exact solution'}, 'Interpreter', 'latex');
    % set(h1,'FontSize',12,'Location','northwest');

end



function da = combofun(t,u) %The inputs u and t are given by the RK4 method
% global lm
A = [0 1 0 0;
    0 0 1 0;
    0 0 0 1;
    -10*u(4) 10*u(3) 10*u(2) -10*(u(1))]; % Jacobian df/dy
dy = [u(2); u(3); u(4); 10*(u(2)*u(3)-u(1)*u(4))]; % This si the value of f=A*y+q. Here y = [u(1); u(2); u(3)]
dz = A*[u(5); u(6); u(7); u(8)]; % This is the value A*zi_j i.e the jth column of zi where j=1,2,3
da = [dy;dz];
end