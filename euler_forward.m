function [t,y] = euler_forward(t0,y0,t_end,n,fcn)
% Solve the initial value problem
% y’ = fcn(t,y), t0 <= t <= b, y(t0)=y0
h = (t_end-t0)/n; % step size.
t = linspace(t0,t_end,n)'; % an array of n points between (and including) t0 and t_end
y = zeros(n,1); % solution vector initialization
y(1) = y0; % initial condition
for i = 2:n
    y(i) = y(i-1)+h*feval(fcn,t(i-1),y(i-1));
end
end

