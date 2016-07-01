%% Produce plots for MVNExample
InitializeWorkspaceDisplay %clean up 
load MVNProbExampleAllData.mat

%% First IID and unscrambled Sobol
axisvec = [100 1e7 1e-12 0.1];
xtick = 10.^(2:7);
ytick = 10.^(-12:2:2);
figure
h = loglog(nvec,errmedMVNProbIIDGn,'.', ...
   [nvec(1) nlarge],errmedMVNProbIIDGn(1)*[1 sqrt(nvec(1)/nlarge)], '--', ...
   [nvec nvec]', [errmedMVNProbIIDGn errtopMVNProbIIDGn]', '-', ...
   'color', MATLABBlue);
hold on
h = [h(1:2); loglog(nvec,errMVNProbuSobolGn,'.', ...
   [nvec(1) nlarge],errMVNProbuSobolGn(1)*[1 nvec(1)/nlarge], '--', ...
   'color', MATLABOrange)];
legend(h,{'IID MC','\(O(n^{-1/2})\)', ...
   'Sobol''','\(O(n^{-1})\)'}, ...
   'location','southwest')
legend boxoff
axis(axisvec)
set(gca,'Xtick',xtick,'YTick',ytick)
xlabel('Sample Size, \(n\)')
ylabel('Error, \(|\mu - \hat{\mu}|\)')
print -depsc MVNIIDUSobol.eps


%% Next, IID, unscrambled Sobol, and scrambled Sobol
figure
h = loglog(nvec,errmedMVNProbIIDGn,'.', ...
   [nvec(1) nlarge],errmedMVNProbIIDGn(1)*[1 sqrt(nvec(1)/nlarge)], '--', ...
   [nvec nvec]', [errmedMVNProbIIDGn errtopMVNProbIIDGn]', '-', ...
   'color', MATLABBlue);
hold on
h = [h(1:2); loglog(nvec,errMVNProbuSobolGn,'.', ...
   [nvec(1) nlarge],errMVNProbuSobolGn(1)*[1 nvec(1)/nlarge], '--', ...
   'color', MATLABOrange)];
h = [h(1:4); loglog(nvec,errmedMVNProbSobolGn,'.', ...
      [nvec(1) nlarge],errmedMVNProbSobolGn(1)*[1 (nvec(1)/nlarge)^1.5], '--', ...
      [nvec nvec]', [errmedMVNProbSobolGn errtopMVNProbSobolGn]' , '-', ...
      'color', MATLABPurple)];
legend(h(1:6),{'IID MC','\(O(n^{-1/2})\)', ...
   'Sobol''','\(O(n^{-1})\)', ...
   'Scrambled Sobol''','\(O(n^{-3/2})\)'}, ...
   'location','southwest')
legend boxoff
axis(axisvec)
set(gca,'Xtick',xtick,'YTick',ytick)
xlabel('Sample Size, \(n\)')
ylabel('Error, \(|\mu - \hat{\mu}|\)')
print -depsc MVNIIDUSobolSobol.eps

%% Plot the Genz function
xplot = (0:0.002:1);
nx = numel(xplot);
[xx,yy] = meshgrid(xplot);
dim = numel(MVNProbIIDGn.a) - 1;
xyplot = [xx(:) yy(:) zeros(nx*nx, dim-2)];
zz = reshape(MVNProbIIDGn.f(xyplot),nx,nx);
figure
surf(xx,yy,zz)
shading interp
xlabel('\(x_1\)')
ylabel('\(x_2\)')
view(-20,10)
print -depsc GenzFun.eps


