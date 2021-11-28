function [error] = FDM_DRICHLET_CONDITIONS(N)
    %To run only this function type into terminal
    %FDM_DRICHLET_CONDITIONS(50) 
    %i.e the name of the function to be run.
    
    % ASR NODES 2020
    % SOLVES y''(t)=exp(x), y(0)=u0, y(1)=u1

    format long
    % N = 50; % no of discritization points
    h = 1/N; % mesh size
    u0=0; % u at left boundary
    u1=0; % u at right boundary
    x = (h:h:1-h)'; % interior nodes
    M = zeros(N-1,N-1);
    b = zeros(N-1,1);

    %  left 
    M(1, 1) = -2/(h^2);
    M(1, 2) = 1/(h^2);
    b(1) = exp(x(1))-u0/(h^2);

    % 
    for i = 2:N-2
        M(i, i-1) = 1/(h^2);
        M(i, i+1) = 1/(h^2);
        M(i, i) = -2/(h^2);
        b(i) = exp(x(i));
    end

    % node on the right
    M(N-1, N-1) = -2/(h^2);
    M(N-1, N-2) = 1/(h^2);
    b(N-1) = exp(x(N-1))-u1/(h^2);

    u = M\b;

    figure(1);plot(x,u,'ok','LineWidth',1),grid on;
    hold on
    plot(x,-1+exp(x)+u0+x-exp(1)*x-u0*x+u1*x,'r','LineWidth',1.5)
    xlabel('$x$','FontSize',13,'Color','k', 'Interpreter', 'latex')
    ylabel('$u(x)$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex')
    error = max(abs(u-(-1+exp(x)+u0+x-exp(1)*x-u0*x+u1*x)));
end