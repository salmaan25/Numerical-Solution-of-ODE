function [y] = classical_RK4()
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
%     title('classical RK4 method, Uniform h');
%     
% hold off
% legend([col1(1),col2(1),col3(1),col4(1),col5(1)], 'n=250', 'n=500', 'n=1000','n=2000', 'n=4000', 'Location','northeast');
% Uncomment these lines for all in one graphs

end

function [t,H] = solveforn(n)
    global M;
    t0 = 0; %Starting point
    t_end = 2000; %Ending point

    y0=[1.4155, 1.545*M/48.888]; %Initial values
    h = (t_end-t0)/n; % step size.
    t = t0:h:t_end; % an array of n points between (and including) t0 and t_end
    y = zeros(length(y0),length(t)); % solution vector initialization length(y0)= 2

    H = zeros(length(t)); % Hamiltonian; should be constant


    y(:,1) = y0; % initial condition
    H0 = feval(@H,y(:,1)); % H(0)

     for i = 2:length(t)
        k_1 = feval(@fcn,t(i-1),y(:,i-1));
        k_2 = feval(@fcn,t(i-1)+0.5*h,y(:,i-1)+0.5*h*k_1);
        k_3 = feval(@fcn,t(i-1)+0.5*h,y(:,i-1)+0.5*h*k_2);
        k_4 = feval(@fcn,t(i-1)+h,y(:,i-1)+h*k_3);

        y(:,i) = y(:,i-1) + (1/6)*h*(k_1+2*k_2+2*k_3+k_4);
        H(i) = feval(@H, y(:,i))-H0;
     end
    
    hold on
        figure
        col1 = plot(t,H, 'Color','m');
        xlabel('t');
        ylabel('H(t)-H(0)');
        title(['classical RK4 method, n = ', num2str(n), ', h = ', num2str(h)]);
    
    hold off
    legend(col1(1), ['n=',num2str(n)],'Location','northeast');
    
    %for plotting y
%     figure
%     plot(t,y);
%     title('p(t) vs t and q(t) vs t');
%     xlabel('t') ;
%     legend({'q', 'p'},'Location','southwest');
end

function [Hpq] = H(y)
global M;
global q0;
global S;
global D;
Hpq = (y(2)^2)/(2*M) + D*( 1-((exp(1))^(-S*(y(1)-q0))) )^2;
end

function [dydt] = fcn(t,y)
global M;
global q0;
global S;
global D;
dydt = [y(2)/M,-2*D*S*( 1-(exp(1))^(-S*(y(1)-q0)) ) * ( (exp(1))^(-S*(y(1)-q0)) )]';
end