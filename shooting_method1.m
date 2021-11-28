%---ASR (NODEs 2020)
% Shooting method Example 7.2 pg 180
% uses classicalRK_4 for IVPs
clc

global lm E Pi B0 Bb
tspan=[0 1]; % set time interval
lm = 2; % lambda
E = 2.718281828459046;
Pi = pi;

B0 = [1, 0, 0;
    0, 0, 0;
    0, 0, 0]; % B0 boundary matrix at t = 0;
Bb = [0, 0, 0;
    1, 0, 0;
    0, 1, 0]; % Bb boundary matrix at t = b;
b1 = 1 + (1 + E^(-2*lm) + E^(-lm))/(2 + E^(-lm));
b2 = 0;
b3 = (3*lm - lm/E^lm)/(2 + E^(-lm));
alpha = [b1 b2 b3]';

I = eye(3);

figure1 = figure('Units','inches','Position',[0,0,6,6*0.7]); 

% first shooting
disp('first shooting')
disp('s is:')
s0  = [b1 0 0]';
disp(s0)

[t,y1]=classicalRK_4(@combofun,tspan, [s0; I(:,1)], 100);
[t,y2]=classicalRK_4(@combofun,tspan, [s0; I(:,2)], 100);
[t,y3]=classicalRK_4(@combofun,tspan, [s0; I(:,3)], 100);


residual = B0*[y1(1,1); y1(2,1); y1(3,1)]+Bb*[y3(1,101); y3(2,101); y3(3,101)]-alpha;
disp('residual is:')
disp(residual)


p0 = plot(t, (E.^(-1 + t)+E.^(2*(-1 + t))+E.^(-t))/(2 + E^(-1))+cos(Pi*t),'k-','LineWidth',1);
hold on

p1 = plot(t,y1(1,:),'g-.');


% second shooting
disp('second shooting')
disp('s is:')

zb = [y1(4:6,101), y2(4:6,101), y3(4:6,101)];
s0 = s0-(B0+Bb*zb)\residual;
disp(s0)

[t,y1]=classicalRK_4(@combofun,tspan, [s0; I(:,1)], 100);
[t,y2]=classicalRK_4(@combofun,tspan, [s0; I(:,2)], 100);
[t,y3]=classicalRK_4(@combofun,tspan, [s0; I(:,3)], 100);


residual = B0*[y1(1,1); y1(2,1); y1(3,1)]+Bb*[y3(1,101); y3(2,101); y3(3,101)]-alpha;
disp('residual is:')
disp(residual)
p2 = plot(t,y1(1,:),'m-.');
hold on

% third shooting
disp('third shooting')
disp('s is:')

zb = [y1(4:6,101), y2(4:6,101), y3(4:6,101)];
s0 = s0-(B0+Bb*zb)\residual;
disp(s0)

[t,y1]=classicalRK_4(@combofun,tspan, [s0; I(:,1)], 100);
[t,y2]=classicalRK_4(@combofun,tspan, [s0; I(:,2)], 100);
[t,y3]=classicalRK_4(@combofun,tspan, [s0; I(:,3)], 100);


residual = B0*[y1(1,1); y1(2,1); y1(3,1)]+Bb*[y3(1,101); y3(2,101); y3(3,101)]-alpha;
disp('residual is:')
disp(residual)

xlabel('$t$','FontSize',13,'Color','k', 'Interpreter', 'latex')
ylabel('$y_1(t)$','FontSize',13,'Color','k', 'Interpreter', 'latex')
title('Shooting method','FontSize',13,'Color','k', 'Interpreter', 'latex')
h1 = legend([p1 p2 p0],{'1st shooting','2nd shooting','exact solution'}, 'Interpreter', 'latex');
set(h1,'FontSize',12,'Location','northwest');





function da = combofun(t,u)
global lm
A = [0 1 0;
    0 0 1;
    -2*lm^2 lm^2 2*lm]; % Jacobian df/dy
dy = A*[u(1); u(2); u(3)]+q(t);
dz = A*[u(4); u(5); u(6)];
da = [dy;dz];
end


function h = q(t)
global lm E Pi 
h = [0, 0, (-2*(E^lm + E^(2*lm*t)+...
    E^(lm*(-1 + 3*t)))*(-1 + lm)*lm^2)/(E^(lm*t)*(1 + 2*E^lm))+... 
  2*lm*(lm + Pi^2)*cos(Pi*t)+Pi*(lm^2 + Pi^2)*sin(Pi*t)]';
end