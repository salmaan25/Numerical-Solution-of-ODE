NN = [10, 20, 40, 80, 160, 320, 640, 1280];
t0 = 0; % initial time a
y0 = 0; % initial condition y(a) = t0
t_end = 1; % final time b
error = zeros(1,length(NN));
for i = 1:length(NN)
    nn = NN(i);
    [t,y_n] = euler_forward(t0,y0,t_end,nn,@fcn);
    y_exact = -t.^2;
    error(i) = max(abs(y_exact-y_n));
end

figure1 = figure('Units','inches','Position',[0,0,6,6*0.7]); 
p1 = loglog(NN,error,'-ok', 'MarkerSize',6,'LineWidth',1);
xlim([10, 1280])
xticks(NN)
xlabel('$N $','FontSize',13,'Color','k', 'Interpreter', 'latex')
ylabel('$\max_{1 \leq n \leq N}|\mathbf{y}(t_n)-\mathbf{y}_n|$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex')
title('Forward Euler (convergence)','FontSize',13,'Color','k', 'Interpreter', 'latex')
saveas(figure1,'C:\Users\aniru\Desktop\recoursematerial\figures\forwardeulerfig2.jpg') 

