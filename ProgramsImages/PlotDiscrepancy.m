%% Discrepancy Plot
gail.InitializeWorkspaceDisplay %clean up 
format long
load discrepancyData

%% Plot geometric interpretation
rng(0)
n=32;
x = net(scramble(sobolset(2),'MatousekAffineOwen'),n);
xcorner = [0.6 0.4];
whichin = (x(:,1) <= 0.6) & (x(:,2) <= 0.4);
figure
plot(x(~whichin,1),x(~whichin,2),'.')
hold on
plot(x(whichin,1),x(whichin,2),'.')
axis square
axis([0 1 0 1])
ticks = 0:0.2:1;
set(gca,'Xtick',ticks,'Ytick',ticks)
plot([0 xcorner([1 1]) 0 0],[0 0 xcorner([2 2]) 0],'k-', ...
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

%% Plot L2 unweighted discrepancy for Sobol' points
compUnWtdDisc = false;
if compUnWtdDisc
   tic
   rng(47)
   mvec = (0:12)';
   nvec = 2.^mvec;
   nmax = max(nvec);
   nn = numel(nvec);
   dvec = [1 4 12 52 250];
   dmax = max(dvec);
   nd = numel(dvec);
   k0 = (4/3).^dvec;
   nrep = 30;
   disc2 = zeros(nn,nd,nrep);

   for ii = 1:nrep
      gail.TakeNote(ii,1)
      x = net(scramble(sobolset(dmax),'MatousekAffineOwen'),nmax);
      Kmat = ones(nmax,nmax);
      kvec = ones(nmax,1);
      kd = 1;
      for k = 1:dmax
         kvec = kvec .* ((3 - x(:,k).^2)/2);
         Kmat = Kmat .* (2 - bsxfun(@max,x(:,k),x(:,k)'));
         if k == dvec(kd)
            temp1 = cumsum(kvec);
            temp2 = cumsum(cumsum(Kmat,1),2);
            disc2(:,kd,ii) = k0(kd) - 2*(temp1(nvec)./nvec) + ...
               temp2(sub2ind([nmax nmax],nvec,nvec))./(nvec.^2);
            kd = kd + 1;
         end
      end
   end
   disc = sqrt(mean(disc2,3));
   scdisc = bsxfun(@rdivide,disc,disc(1,:));
   toc
end
figure
h = loglog(nvec,scdisc,'.');
hold on
h = [h; loglog(nvec([1 nn]),scdisc(1,1)*[1 nvec(1)/nvec(nn)],'--','color',get(h(1),'color'))];
%h = [h; loglog(nvec([1 nn]),scdisc(1,nd)*[1 sqrt(nvec(1)/nvec(nn))],'--','color',get(h(nd),'color'))];
h = [h; loglog(nvec([1 nn]),[1 sqrt(nvec(1)/nvec(nn))],'k-')];
xlabel('\(n\)')
ylabel('Scaled DSC')
legendstuff = cell(1,nd+2);
legendstuff{1} = ['\(d = ' int2str(dvec(1)) '\)'];
legendstuff{2} = '\(O(n^{-1})\)';
for k = 2:nd
   legendstuff{k+1} = ['\(d = ' int2str(dvec(k)) '\)'];
end
%legendstuff{nd+2} = '\(O(n^{-1/2})\)';
legendstuff{nd+2} = 'IID, \(O(n^{-1/2})\)';
legend(h([1 nd+1 2:nd nd+2]),legendstuff,'location','southwest')
legend boxoff
title('Sobol'' Points')
print -depsc UnwtL2Disc.eps

%% Plot L2 weighted discrepancy for Sobol' points
CompWtdDisc = false;
if CompWtdDisc
   tic
   rng(47)
   gamma = (1:dmax).^-3;
   k0 = cumprod(1 + gamma/3);
   k0 = k0(dvec);
   nrep = 30;
   wtddisc2 = zeros(nn,nd,nrep);

   for ii = 1:nrep
      gail.TakeNote(ii,1)
      x = net(scramble(sobolset(dmax),'MatousekAffineOwen'),nmax);
      Kmat = ones(nmax,nmax);
      kvec = ones(nmax,1);
      kd = 1;
      for k = 1:dmax
         kvec = kvec .* ((1 + gamma(k)/2) - (gamma(k)/2) * x(:,k).^2);
         Kmat = Kmat .* ((1+gamma(k)) - gamma(k)*bsxfun(@max,x(:,k),x(:,k)'));
         if k == dvec(kd)
            temp1 = cumsum(kvec);
            temp2 = cumsum(cumsum(Kmat,1),2);
            wtddisc2(:,kd,ii) = k0(kd) - 2*(temp1(nvec)./nvec) + ...
               temp2(sub2ind([nmax nmax],nvec,nvec))./(nvec.^2);
            kd = kd + 1;
         end
      end
   end
   wtddisc = sqrt(mean(wtddisc2,3));
   scwtddisc = bsxfun(@rdivide,wtddisc,wtddisc(1,:));
   toc
end
figure
h = loglog(nvec,scwtddisc,'.');
hold on
h = [h; loglog(nvec([1 nn]),scwtddisc(1,1)*[1 nvec(1)/nvec(nn)],'--','color',get(h(1),'color'))];
xlabel('\(n\)')
ylabel('Scaled, Weighted DSC')
legendstuff = cell(1,nd+1);
legendstuff{1} = ['\(d = ' int2str(dvec(1)) '\)'];
legendstuff{2} = '\(O(n^{-1})\)';
for k = 2:nd
   legendstuff{k+1} = ['\(d = ' int2str(dvec(k)) '\)'];
end
legend(h([1 nd+1 2:nd]),legendstuff,'location','southwest')
legend boxoff
title('Sobol'' Points')
print -depsc WtL2Disc.eps

%% Plot L2 unweighted discrepancy for lattice points
compUnWtdDiscLat = true;
if compUnWtdDiscLat
   tic
   rng(47)
   mvec = (0:12)';
   nvec = 2.^mvec;
   nmax = max(nvec);
   nn = numel(nvec);
   dvec = [1 4 12 52 250];
   dmax = max(dvec);
   nd = numel(dvec);
   k0 = (4/3).^dvec;
   nrep = 30;
   disc2 = zeros(nn,nd,nrep);

   for ii = 1:nrep
      gail.TakeNote(ii,1)
      x = mod(bsxfun(@plus,gail.lattice_gen(1,nmax,dmax),rand(1,dmax)),1);
      Kmat = ones(nmax,nmax);
      kvec = ones(nmax,1);
      kd = 1;
      for k = 1:dmax
         kvec = kvec .* ((3 - x(:,k).^2)/2);
         Kmat = Kmat .* (2 - bsxfun(@max,x(:,k),x(:,k)'));
         if k == dvec(kd)
            temp1 = cumsum(kvec);
            temp2 = cumsum(cumsum(Kmat,1),2);
            disc2(:,kd,ii) = k0(kd) - 2*(temp1(nvec)./nvec) + ...
               temp2(sub2ind([nmax nmax],nvec,nvec))./(nvec.^2);
            kd = kd + 1;
         end
      end
   end
   disc = sqrt(mean(disc2,3));
   scdisclat = bsxfun(@rdivide,disc,disc(1,:));
   toc
end
figure
h = loglog(nvec,scdisclat,'.');
hold on
h = [h; loglog(nvec([1 nn]),scdisclat(1,1)*[1 nvec(1)/nvec(nn)],'--','color',get(h(1),'color'))];
%h = [h; loglog(nvec([1 nn]),scdisc(1,nd)*[1 sqrt(nvec(1)/nvec(nn))],'--','color',get(h(nd),'color'))];
h = [h; loglog(nvec([1 nn]),[1 sqrt(nvec(1)/nvec(nn))],'k-')];
xlabel('\(n\)')
ylabel('Scaled DSC')
legendstuff = cell(1,nd+2);
legendstuff{1} = ['\(d = ' int2str(dvec(1)) '\)'];
legendstuff{2} = '\(O(n^{-1})\)';
for k = 2:nd
   legendstuff{k+1} = ['\(d = ' int2str(dvec(k)) '\)'];
end
%legendstuff{nd+2} = '\(O(n^{-1/2})\)';
legendstuff{nd+2} = 'IID, \(O(n^{-1/2})\)';
legend(h([1 nd+1 2:nd nd+2]),legendstuff,'location','southwest')
legend boxoff
title('Lattice Points')
print -depsc UnwtL2Disclat.eps

%% Plot L2 weighted discrepancy for lattice points
CompWtdDiscLat = true;
if CompWtdDiscLat
   tic
   rng(47)
   gamma = (1:dmax).^-3;
   k0 = cumprod(1 + gamma/3);
   k0 = k0(dvec);
   nrep = 30;
   wtddisc2 = zeros(nn,nd,nrep);

   for ii = 1:nrep
      gail.TakeNote(ii,1)
      x = mod(bsxfun(@plus,gail.lattice_gen(1,nmax,dmax),rand(1,dmax)),1);
      Kmat = ones(nmax,nmax);
      kvec = ones(nmax,1);
      kd = 1;
      for k = 1:dmax
         kvec = kvec .* ((1 + gamma(k)/2) - (gamma(k)/2) * x(:,k).^2);
         Kmat = Kmat .* ((1+gamma(k)) - gamma(k)*bsxfun(@max,x(:,k),x(:,k)'));
         if k == dvec(kd)
            temp1 = cumsum(kvec);
            temp2 = cumsum(cumsum(Kmat,1),2);
            wtddisc2(:,kd,ii) = k0(kd) - 2*(temp1(nvec)./nvec) + ...
               temp2(sub2ind([nmax nmax],nvec,nvec))./(nvec.^2);
            kd = kd + 1;
         end
      end
   end
   wtddisc = sqrt(mean(wtddisc2,3));
   scwtddisclat = bsxfun(@rdivide,wtddisc,wtddisc(1,:));
   toc
end
figure
h = loglog(nvec,scwtddisclat,'.');
hold on
h = [h; loglog(nvec([1 nn]),scwtddisclat(1,1)*[1 nvec(1)/nvec(nn)],'--','color',get(h(1),'color'))];
xlabel('\(n\)')
ylabel('Scaled, Weighted DSC')
legendstuff = cell(1,nd+1);
legendstuff{1} = ['\(d = ' int2str(dvec(1)) '\)'];
legendstuff{2} = '\(O(n^{-1})\)';
for k = 2:nd
   legendstuff{k+1} = ['\(d = ' int2str(dvec(k)) '\)'];
end
legend(h([1 nd+1 2:nd]),legendstuff,'location','southwest')
legend boxoff
title('Lattice Points')
print -depsc WtL2Disclat.eps

%% Save Data
save discrepancyData

