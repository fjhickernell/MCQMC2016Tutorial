%% Test out MLE
function [muhat,aMLE,err,out] = TestOutMLE(dim,BernPolynOrder,ptransform)

%{
if nargin < 3
   f = @(x) cos(x(:,1) + x(:,2)); exactInteg = -cos(2) + 2*cos(1) - 1;
   if nargin < 2
      dim = 2;
      if nargin < 1
         nvec = 2.^(10:20); %10:20
      end
   end
end
%}

nvec = 2.^(8:20); %20

domain = repmat([0;1],[1,dim]);

if dim==2
    f = @(x) cos(x(:,1) + x(:,2)); exactInteg = 4 *cos(1) * sin(1/2)^2;
    
    %f = @(x) sin(pi*x(:,1)).*sin(pi*x(:,2)); exactInteg = 4/pi^2;
    
    %f = @(x) cos(2*pi*x(:,1).^2).*cos(2*pi*x(:,2).^2); exactInteg = 0.0595978471360429;
    fName='Exp(Cos)'; ffunc = @(x,y) exp(cos(x)+cos(y)); f = @(x) ffunc(x(:,1),x(:,2)); exactInteg = 5.48297273934351;
end
if dim==3
    %coeff = [1 (1/2)^4 (1/2)^8];
    coeff = [1 1 1];
    %coeff = [2*pi 2*pi 2*pi];
    a=coeff(1);b=coeff(2);c=coeff(3);
    
    fName='Cos'; f = @(x) cos(a*x(:,1) + b*x(:,2) + c*x(:,3)); exactInteg = (8 *cos((a+b+c)/2) *sin(a/2) *sin(b/2)*sin(c/2))/(a*b*c); 
    %f = @(x) cos(a*x(:,1).^2 + b*x(:,2).^2 + c*x(:,3).^2); 
    
    %f = @(x) cos(2*pi*x(:,1).^2).* cos(2*pi*x(:,2).^2) .* cos(2*pi*x(:,3).^2); exactInteg = 0.0145494259294652;
    %ffunc = @(x,y,z) exp(cos(x)+cos(y)+cos(z)); f = @(x) ffunc(2*pi*x(:,1),2*pi*x(:,2),2*pi*x(:,3)); exactInteg = 2.02940587037004;
    %fName='Exp(Cos)'; ffunc = @(x,y,z) exp(cos(x)+cos(y)+cos(z)); f = @(x) ffunc(a*x(:,1),b*x(:,2),c*x(:,3)); exactInteg = 12.8387910242451;
    
    %exactInteg =  8 *cos(3/2) * sin(1/2)^3;
    %
    
    fName='Exp(Cos)'; ffunc = @(x,y,z) exp(cos(2*pi*x)+cos(2*pi*y)+cos(2*pi*z)); f = @(x) ffunc(x(:,1),x(:,2),x(:,3)); exactInteg = 2.0294058703700;
    
    
    %f = @(x) sin(pi*x(:,1)).*sin(pi*x(:,2)).*sin(pi*x(:,3)); exactInteg = 8/pi^3;
end
if dim==4
    %coeff = [1 (1/2)^4 (1/2)^8 (1/2)^12];
    coeff = [1 1 1 1];
    a=coeff(1);b=coeff(2);c=coeff(3);d=coeff(4);
    fName='Cos'; f = @(x) cos(a*x(:,1) + b*x(:,2) + c*x(:,3) + d*x(:,4)); 
    %exactInteg =  16 *cos(2) * sin(1/2)^4;
    
    exactInteg = 8*sin(b/2)*sin(c/2)*sin(d/2)*...
        (sin((2*a+b+c+d)/2) - sin((b+c+d)/2))/(a*b*c*d);

end

% General form
% f = @(x) cos(sum(x,2)); exactInteg =  2^dim *cos(dim/2) * sin(1/2)^dim;
%


%fName='\(\cos(\sum_{i=1}^d x_i)\)';
figSavePath = '/home/jagadees/MyWriteup/Mar1stweek/';
%cubMLE(f,nvec,domain,whSample,whKer,powerFuncMethod)
whSample = 'Lattice1';
whKer = 'Mat1';
whKer = 'Fourier';
powerFuncMethod = 'Cauchy';

fullPath = strcat(figSavePath,'/',fName,'/',ptransform,'/');
if exist(fullPath,'dir')==false
    mkdir(fullPath);
end

[muhat,out] = cubMLE(f,nvec,domain, whSample,whKer,powerFuncMethod,BernPolynOrder,ptransform,fName,fullPath);
aMLE = out.aMLE
err = abs(exactInteg - muhat)
errbd = out.ErrBd
disc = out.disc2
time = out.time
%BernPolynOrder = out.BernPolynOrder;

%% plot error
errCubMLE = abs(exactInteg - muhat)

hFigErr = figure;MATLABYellow = [0.9290, 0.6940, 0.1250];
set(hFigErr, 'units', 'inches', 'Position', [4 4 10 7])
loglog(nvec,errCubMLE,'.', ...
    nvec,abs(out.ErrBd),'-.', ...
      [nvec(1) nvec(end)],errCubMLE(1)*[1 (nvec(1)/nvec(end))^2], '--', ...
      [nvec(1) nvec(end)],errCubMLE(1)*[1 (nvec(1)/nvec(end))^1], ':', ...
      'color', MATLABYellow);
legend({'Bayesian Cubature Fourier Lattice Actual Error', ...
    'Bayesian Cubature Fourier Lattice Error Bound', '\(O(n^{-2})\)', '\(O(n^{-1})\)'}, ...
   'location','southwest')
legend boxoff
xlabel('Sample Size, \(n\)')
ylabel('Error, \(|\mu - \hat{\mu}|\)')
axis tight

title(sprintf('%s d=%d Bernoulli=%d, PeriodTx=%s', fName, dim, BernPolynOrder, out.ptransform))
saveas(hFigErr, sprintf('%s%s Error d_%d bernoulli_%d Period_%s.png', fullPath, fName, dim, BernPolynOrder, out.ptransform))
fprintf('done');


end