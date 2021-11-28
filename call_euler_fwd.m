t0 = 0; % initial time a
y0 = 0; % initial condition y(a) = t0
t_end = 1; % final time b
[t1,y1] = euler_forward(t0,y0,t_end,100,@fcn); %n=100; % number of points N
[t2,y2] = euler_forward(t0,y0,t_end,10,@fcn); %n=10; % number of points N
%Anylatic solution
y1_exact = -t1.^2;
y2_exact = -t2.^2;
%results: post-processing 
figure1 = figure('Units','inches','Position',[0,0,6,6*0.7]); 
p1 = plot(t1,y1,'ok',t2,y2,'sm', 'MarkerSize',6,'LineWidth',1);
% set(p1,'Color',[0 0 0])
hold on
p2 = plot(t1,y1_exact,'LineWidth',1);
set(p2,'Color',[1 0 0])
xlabel('$t$','FontSize',13,'Color','k', 'Interpreter', 'latex')
ylabel('$\mathbf{y}$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex')
title('Forward Euler','FontSize',13,'Color','k', 'Interpreter', 'latex')
h1 = legend([p1(1) p1(2) p2],{'N = 100','N=10','Exact Solution'}, 'Interpreter', 'latex');
set(h1,'FontSize',12,'Location','southwest');
%print -depsc 'C:\Users\aniru\Desktop\recoursematerial\figures\forwardeulerfig1.jpg'
saveas(figure1,'C:\Users\aniru\Desktop\recoursematerial\figures\forwardeulerfig1.jpg') 
% Error plot
figure2 = figure('Units','inches','Position',[0,0,6,6*0.7]);
p1 = semilogy(t1,abs(y1-y1_exact),'ok',t2,abs(y2-y2_exact),'sm', 'MarkerSize',6,'LineWidth',1);
xlabel('$t$','FontSize',13,'Color','k', 'Interpreter', 'latex')
ylabel('$|\mathbf{y}(t_n)-\mathbf{y}_n|$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex')
title('Forward Euler (error plot)','FontSize',13,'Color','k', 'Interpreter', 'latex')
h1 = legend([p1(1) p1(2)],{'N = 100','N=10'}, 'Interpreter', 'latex');
set(h1,'FontSize',12,'Location','southeast');
%print -depsc 'C:\Users\aniru\Desktop\recoursematerial\figures\forwardeulerfig2.jpg'
saveas(figure2,'C:\Users\aniru\Desktop\recoursematerial\figures\forwardeulerfig2.jpg') 
