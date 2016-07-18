%% Discrepancy Plot
InitializeWorkspaceDisplay %clean up 
format long

%% Plot geometric interpretation
rng(0)
n=32;
x = net(scramble(sobolset(2),'MatousekAffineOwen'),n);
xcorner = [0.6 0.4];
figure
plot(x(:,1),x(:,2),'.')
hold on
plot([xcorner([1 1]) 0],[0 xcorner([2 2])],'k-', ...
   xcorner(1), xcorner(2), '.k')
xlabel('\(x_5\)')
ylabel('\(x_8\)')
text(xcorner(1), xcorner(2)+0.05, ...
   ['\((' num2str(xcorner(1)) ', '  num2str(xcorner(2)) ')\)'])
localdisc = prod(xcorner) ...
   - sum((x(:,1) <= xcorner(1)) & (x(:,2) <= xcorner(2)))/n
print -depsc LocalDiscrep.eps

%% Plot kernelfigure
x = (0:0.002:1)';
[xx, yy] = meshgrid(x);
Kxy = 2 - abs(xx-yy);
figure
surf(xx,yy,Kxy)
shading interp
hold on 
t = 0.3;
h = plot3(x,t*ones(size(x)),2 - abs(x-t),'-k');
xlabel('\(x\)')
ylabel('\(t\)')
zlabel('\(K(x,t)\)')
xtick = 0:0.2:1;
ztick = 1:0.2:2;
%set(gca,'Xtick',xtick,'Ytick',xtick,'Ztick',ztick)
view(-20,20)
legend(h,['\(K(x,' num2str(t) ')\)'],'Location','northwest')
legend boxoff
print -depsc L2Kernel.eps

