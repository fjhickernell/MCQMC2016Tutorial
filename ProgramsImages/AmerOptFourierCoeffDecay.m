function AmerOptFourierCoeffDecay(testfun,d,printc)
%[~,~,~,MATLABVERSION] = GAILstart(false);

%if usejava('jvm') || MATLABVERSION <= 7.12
    %% Garbage collection and initialization
    format compact %remove blank lines from output
    format long e %lots of digits
    %clearvars %clear all variables
    close all %close all figures
    ColorOrder=get(gca,'ColorOrder'); close all
    set(0,'defaultaxesfontsize',24,'defaulttextfontsize',24) %make font larger
    set(0,'defaultLineLineWidth',3) %thick lines
    set(0,'defaultTextInterpreter','latex') %latex axis labels
    set(0,'defaultLegendInterpreter','latex') %latex axis labels
    set(0,'defaultLineMarkerSize',40) %latex axis labels
    if nargin < 3; printc='color'; end

    %% Initialize parameters
    mmax=22; %maximum number of points is 2^mmax
    mdualvec=12:12;
    mplot=22;
    mlag=4;
    %testfun=@(x) exp(-3*x).*sin(10*x.^2); d=1; %test function
    sobstr=sobolset(d);
    sobstr=scramble(sobstr,'MatousekAffineOwen');

%     %% Plot function
%     figure
%     xplot=(0:0.002:1);
%     yplot=testfun(xplot);
%     plot(xplot,yplot,'-');
%     yminf=1.1*min(yplot);
%     ymaxf=1.1*max(yplot);
%     axis([0 1 yminf ymaxf])
%     xlabel('$x$','interpreter','latex')
%     ylabel('$f(x)$','interpreter','latex')
%     print('FunctionWalshFourierCoeffDecay','-depsc');

    %% Evaluate Function and FWT
    n=2^mmax;
    xpts=sobstr(1:n,1:d);
    y=testfun(xpts);
    %yval=y;
    %yfwt=fwht(y);

    %% Compute FWT
    for l=0:mmax-1
       nl=2^l;
       nmmaxlm1=2^(mmax-l-1);
       ptind=repmat([true(nl,1); false(nl,1)],nmmaxlm1,1);
       evenval=y(ptind);
       oddval=y(~ptind);
       y(ptind)=(evenval+oddval)/2;
       y(~ptind)=(evenval-oddval)/2;
    end

    %% Create kappanumap
    kappanumap=(1:n)'; %initialize map
    for l=mmax-1:-1:1
      nl=2^l;
      oldone=abs(y(kappanumap(2:nl))); %earlier values of kappa, don't touch first one
      newone=abs(y(kappanumap(nl+2:2*nl))); %later values of kappa, 
      flip=find(newone>oldone); %which in the pair are the larger ones
      if ~isempty(flip)
          flipall=bsxfun(@plus,flip,0:2^(l+1):2^mmax-1);
          flipall=flipall(:);
          temp=kappanumap(nl+1+flipall); %then flip 
          kappanumap(nl+1+flipall)=kappanumap(1+flipall); %them
          kappanumap(1+flipall)=temp; %around
      end
    end
    ymap=y(kappanumap);  

    %% Plot FW coefficients
    ltgray=0.8*ones(1,3);
    gray=0.5*ones(1,3);
    nplot=2^mplot;
    yfwtabs=abs(ymap(1:nplot));
    ymin=max(1e-15,min(yfwtabs));
    ymax=max([1; yfwtabs]);
    for mdual=mdualvec
       ndual=2^mdual;
       whdual=ndual*(1:2^(mplot-mdual)-1);
       whsmall=1:ndual-1;
       whbig=ndual:nplot-1;
       muse=mdual-mlag;
       nuse=2^muse;
       whuse=nuse/2:nuse-1;
       figure
       errbdsum = sum(abs(yfwtabs(whuse+1)))
       errdualsum = sum(abs(yfwtabs(whdual+1)))
       inflate = 5/nuse
       errbdguess = inflate * errbdsum
       switch printc
           case 'color'
               h=loglog(whsmall,yfwtabs(whsmall+1),'.',...
                  whbig,yfwtabs(whbig+1),'.',...
                  whuse,yfwtabs(whuse+1),'.',...
                  whdual,yfwtabs(whdual+1),'.','MarkerSize',10);
               set(h([3 4]),'MarkerSize',20)
               set(h(1),'Color','k');
               set(h(2),'Color',0.7*[1 1 1]);
               set(h(3),'Color',ColorOrder(5,:));
               set(h(4),'Color',ColorOrder(2,:));
           case 'bw'
               h=zeros(4,1);
               h(1)=loglog(whsmall,yfwtabs(whsmall+1),'.',...
                  'MarkerSize',10,'MarkerFaceColor',ltgray,'MarkerEdgeColor',ltgray);
               hold on
               h(2)=loglog(whbig,yfwtabs(whbig+1),'.','MarkerSize',10,...
                  'MarkerFaceColor',gray,'MarkerEdgeColor',gray);
               h(3)=loglog(whuse,yfwtabs(whuse+1),'sk','MarkerSize',7,...
                  'MarkerFaceColor','k');
               h(4)=loglog(whdual,yfwtabs(whdual+1),'.k','MarkerSize',30);
       end
       maxexp=floor(log10(nplot-1));
       set(gca,'Xtick',10.^(0:maxexp))
       axis([1 nplot-1 ymin ymax])
       xlabel('$\kappa$','interpreter','latex')
       ylabel('$|\hat{f}_{\kappa}|$','interpreter','latex')
       legend(h([4 3]),{'error bound',...
          ['$S_{' int2str(mdual) '}(f)$']},...
          'location','southwest','interpreter','latex')
       legend('boxoff')
       set(gca,'Position',[0.2 0.155 0.75 0.77])
       print(['WalshFourierCoeffDecay' int2str(nuse) '.eps'], '-depsc');
    end

    return
    %Plot filtered function bad function
    ydwt=y; %save original coefficients
    whkill=kappanumap(2^7+1:2^8); %choose ones to shrink
    y(whkill)=1e-6*y(whkill); %shrink them
    yfwtabs=abs(y(kappanumap));
    figure
    loglog(1:nplot-1,yfwtabs(2:nplot),'k.','MarkerSize',10);
    maxexp=floor(log10(nplot-1));
    set(gca,'Xtick',10.^(0:maxexp))
    axis([1 nplot-1 ymin ymax])
    xlabel('$\kappa$','interpreter','latex')
    ylabel('$|\hat{f}_{\kappa}|$','interpreter','latex')
    set(gca,'Position',[0.2 0.155 0.75 0.77])
    print('WalshFourierCoeffDecayFilter.eps', '-depsc');
    
    %% Compute FWT
    for l=0:mmax-1
       nl=2^l;
       nmmaxlm1=2^(mmax-l-1);
       ptind=repmat([true(nl,1); false(nl,1)],nmmaxlm1,1);
       evenval=y(ptind);
       oddval=y(~ptind);
       y(ptind)=(evenval+oddval)/2;
       y(~ptind)=(evenval-oddval)/2;
    end
    yfilter=(2.^mmax)*y;
    nplotf=2^10;
    [xsort,ii]=sort(xpts(1:nplotf));
    figure; 
    h=plot(xsort,yfilter(ii(1:nplotf)),'-');
    set(h,'Color',ColorOrder(2,:));
    xlabel('$x$','interpreter','latex')
    ylabel('$f(x)$','interpreter','latex')
    axis([0 1 yminf ymaxf])
    print('FilteredFunctionWalshFourierCoeffDecay', '-depsc');
    err=abs(testfun(xsort)-yfilter(ii(1:nplotf)));
    figure; plot(xsort,err)

    %close all
end
%end



