%% Discrepancy Plot
InitializeWorkspaceDisplay %clean up 
format long
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
