%% Keister's Example of Multidimensional Integration
%
% B. D. Keister, Multidimensional quadrature algorithms, _Computers in
% Physics_, *10*, pp. 119-122, 1996, presents the following
% multidimensional integral, inspired by a physics application:
%
% \[ I = \int_{\mathbb{R}^d} \cos(\lVert \boldsymbol{x} \rVert)
% \exp(-\lVert \boldsymbol{x} \rVert^2) \, \mathrm{d} \boldsymbol{x},
% \qquad d = 1, 2, \ldots. \]

function [muhat,aMLE,err,out] = KeisterCubatureExampleBayesian(dim,BernPolyOrder,ptransform)

normsqd = @(t) sum(t.*t,2); %squared l_2 norm of t

domain = repmat([0;1],[1,dim]);
nvec = 2.^(10:20);

replaceZeros = @(t) (t+(t==0)*eps); % to avoid getting infinity, NaN
yinv = @(t)(erfcinv( replaceZeros(abs(t)) ));  %using erfcinv is more accurate than erfinv with -1
fKeister = @(t,dim) cos( sqrt( normsqd(yinv(t)) )) *(sqrt(pi)/2)^dim;


%% Bayesian Cubature
fName='Keister';
figSavePath = '/home/jagadees/MyWriteup/Apr1stweek/';
whSample = 'Lattice1';
%whKer = 'Mat1';
whKer = 'Fourier';
powerFuncMethod = 'Cauchy';
f1 = @(x) fKeister(x,dim);

fullPath = strcat(figSavePath,'/',fName,'/',ptransform,'/');
if exist(fullPath,'dir')==false
    mkdir(fullPath);
end

[muhat,out] = cubMLE(f1,nvec,domain, whSample,whKer,powerFuncMethod,BernPolyOrder,ptransform,fName,fullPath);


%% plot error
exactInteg = Keistertrue(dim);
errCubMLE = abs(exactInteg - muhat)

plotCubatureError(dim, nvec, errCubMLE, out.ErrBd, fName, BernPolyOrder, out.ptransform, fullPath)


fprintf('Done')
