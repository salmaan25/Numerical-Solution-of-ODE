function [] = RKDP_45()
% y(2) = p, y(1) = q;
% Initializing constants
global M;
M = 0.9953;
global q0;
q0 = 1.41;
global S;
S = 1.814;
global D;
D = 0.0378652;

[t1,H1] = solveforn(250);
[t2,H2] = solveforn(500);
[t3,H3] = solveforn(1000);
[t4,H4] = solveforn(2000);
[t5,H5] = solveforn(4000);

% Uncomment these lines for all in one graphs
% hold on
%     col1 = plot(t1,H1, 'Color','r');
%     col2 = plot(t2,H2, 'Color','g');
%     col3 = plot(t3,H3, 'Color','b');
%     col4 = plot(t4,H4, 'Color','c');
%     col5 = plot(t5,H5, 'Color','m');
%     xlabel('t');
%     ylabel('H(t)-H(0)');
%     title('RKDP45 method, Uniform h');
%     
% hold off
% legend([col1(1),col2(1),col3(1),col4(1),col5(1)], 'n=250', 'n=500', 'n=1000','n=2000', 'n=4000', 'Location','northeast');
% Uncomment these lines for all in one graphs

end

function [t,H] = solveforn(n)
    global M;
    t0 = 0; %Starting point
    t_end = 2000; %Ending point

    y0=[1.4155, 1.545*M/48.888];
    h = (t_end-t0)/n; % step size.
    t = t0:h:t_end; % an array of n points between (and including) t0 and t_end
    y = zeros(length(y0),length(t)); % solution vector initialization length(y0)= 2
    ys = zeros(length(y0),length(t)); % solution vector initialization

    H = zeros(length(t)); % Hamiltonian; should be constant

    y(:,1) = y0; % initial condition
    ys(:,1) = y0; % initial condition
    H0 = feval(@H,y(:,1));

    % J. R. Dormand and P. J. Prince coeffcients
    A = [0, 0, 0, 0, 0, 0, 0;
        1/5, 0, 0, 0, 0, 0, 0;
        3/40, 9/40, 0, 0, 0, 0, 0;
        44/45, -56/15, 32/9, 0, 0, 0, 0;
        19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0, 0;
        9017/3168, -355/33, 46732/5247, 49/176 , -5103/18656, 0, 0;
        35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0];


    b1 = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40]';
    b2 = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0]';

    c = [0, 1/5, 3/10, 4/5, 8/9, 1, 1]';


    for i = 2:length(t)
        t(i)=t(i-1)+h;

        k_1 = feval(@fcn,t(i-1), y(:,i-1));
        k_2 = feval(@fcn,t(i-1)+c(2)*h,y(:,i-1)+h*A(2,1)*k_1);
        k_3 = feval(@fcn,t(i-1)+c(3)*h,y(:,i-1)+h*A(3,1)*k_1+h*A(3,2)*k_2);
        k_4 = feval(@fcn,t(i-1)+c(4)*h,y(:,i-1)+h*A(4,1)*k_1+h*A(4,2)*k_2+h*A(4,3)*k_3);
        k_5 = feval(@fcn,t(i-1)+c(5)*h,y(:,i-1)+h*A(5,1)*k_1+h*A(5,2)*k_2+h*A(5,3)*k_3+h*A(5,4)*k_4);
        k_6 = feval(@fcn,t(i-1)+c(6)*h,y(:,i-1)+h*A(6,1)*k_1+h*A(6,2)*k_2+h*A(6,3)*k_3+h*A(6,4)*k_4+h*A(6,5)*k_5);
        k_7 = feval(@fcn,t(i-1)+c(7)*h,y(:,i-1)+h*A(7,1)*k_1+h*A(7,2)*k_2+h*A(7,3)*k_3+h*A(7,4)*k_4+h*A(7,5)*k_5+h*A(7,6)*k_6);

        y(:,i) = y(:,i-1)+h*b1(1)*k_1+h*b1(2)*k_2+h*b1(3)*k_3+h*b1(4)*k_4+h*b1(5)*k_5+h*b1(6)*k_6+h*b1(7)*k_7;
        ys(:,i) = ys(:,i-1)+h*b2(1)*k_1+h*b2(2)*k_2+h*b2(3)*k_3+h*b2(4)*k_4+h*b2(5)*k_5+h*b2(6)*k_6+h*b2(7)*k_7;
        H(i) = feval(@H, y(:,i))-H0;

    end
    
    hold on
        figure
        col1 = plot(t,H, 'Color','m');
        xlabel('t');
        ylabel('H(t)-H(0)');
        title(['RKDP45 method, n = ', num2str(n), ', h = ', num2str(h)]);
    
    hold off
    legend(col1(1), ['n=',num2str(n)],'Location','northeast');
    
    %for plotting y
%     figure
%     plot(t,y);
%     title('p(t) vs t and q(t) vs t');
%     xlabel('t') ;
%     legend({'q', 'p'},'Location','southwest');
end

function [] = checkWithode45(t0,t_end,y0,H0)
    [tc,yc] = ode45(@fcn,[t0,t_end],y0);
%     fprintf("%d %d\n", length(tc), length(yc(1)));
    HC = zeros(length(tc));
    for i = 1:length(tc)
%         fprintf("%d\n",i);
        HC(i) = feval(@H, yc(1,:))-H0;
    end
    figure
    plot(tc,HC);
end
function [Hpq] = H(y)
% Initializing constants
global M;
global q0;
global S;
global D;
Hpq = (y(2)^2)/(2*M) + D*( 1-((exp(1))^(-S*(y(1)-q0))) )^2;
end

function [dydt] = fcn(t,y)
% Initializing constants
global M;
global q0;
global S;
global D;
% fprintf("eval values p0 = %f, q0 = %f\n", y(2), y(1));
dydt = [y(2)/M,-2*D*S*( 1-(exp(1))^(-S*(y(1)-q0)) ) * ( (exp(1))^(-S*(y(1)-q0)) )]';
end