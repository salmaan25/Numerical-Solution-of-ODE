function [] = sympletic_euler()
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
[t1,H1] = solveforh(2);
[t2,H2] = solveforh(2.3684);
[t3,H3] = solveforh(2.3685);

% Uncomment these lines for all in one graphs
% hold on
%     col1 = plot(t1,H1, 'Color','r');
%     col2 = plot(t2,H2, 'Color','g');
%     col3 = plot(t3,H3, 'Color','b');
%     xlabel('t');
%     ylabel('H(t)-H(0)');
%     title('Sympletic Euler method N =1000');
%     
% hold off
% legend([col1(1),col2(1),col3(1)], 'h=2', 'h=2.3684', 'h=2.3685','Location','northeast');
% Uncomment these lines for all in one graphs

end

function [t,H] = solveforh(h)
    global M;
    t0 = 0; %Starting point
    n = 1000;
    t_end = n*h; %Ending point

    y0=[1.4155, 1.545*M/48.888]; %Initial values
%     h = 2.3684; % step size
%     n = round((t_end-t0)/h); %Number of iterations
    t = t0:h:t_end; % an array of n points between (and including) t0 and t_end
    y = zeros(length(y0),length(t)); % solution vector initialization
    H = zeros(length(t)); % Hamiltonian-H(0); should be 0


    y(:,1) = y0; % initial condition
    H0 = feval(@H,y(:,1)); % H(0)

    for i = 2:n
        y(:,i) = y(:,i-1)+h*feval(@fcn,t(i-1),y(:,i-1),h); %(q_i,p_i)
        H(i) = feval(@H, y(:,i))-H0; %H(i)-H0
    end

    hold on
        figure
        col1 = plot(t,H, 'Color','r');
        xlabel('t');
        ylabel('H(t)-H(0)');
        title(['Sympletic Euler method, N = 1000, h = ', num2str(h)]);
    
    hold off
    legend(col1(1), ['h=',num2str(h)],'Location','northeast');
    
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

function [dydt] = fcn(t,y,h)
global M;
global q0;
global S;
global D;
v1 = (1/M)*(y(2)+h*(-2*D*S*( 1-(exp(1))^(-S*(y(1)-q0)) ) * ( (exp(1))^(-S*(y(1)-q0)) )));
dydt = [v1,-2*D*S*( 1-(exp(1))^(-S*(y(1)-q0)) ) * ( (exp(1))^(-S*(y(1)-q0)) )]';
end