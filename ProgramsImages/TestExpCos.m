%% Test out MLE
function [nvec,muhat,aMLE,errCubMLE,out] = TestExpCos(dim,BernPolyOrder,ptransform)

nvec = 2.^(8:20);

domain = repmat([0;1],[1,dim]);

f = @(x) exp(sum(cos(2*pi*x), 2));
exactInteg = besseli(0,1)^dim;
fName = 'Exp(cos)';
figSavePath = '/home/jagadees/MyWriteup/Apr1stweek/';

whSample = 'Lattice1';
whKer = 'Fourier';
powerFuncMethod = 'Cauchy';

fullPath = strcat(figSavePath,'/',fName,'/',ptransform,'/');
if exist(fullPath,'dir')==false
    mkdir(fullPath);
end

[muhat,out] = cubMLE(f,nvec,domain, whSample,whKer,powerFuncMethod,BernPolyOrder,...
    ptransform,fName,fullPath);
aMLE = out.aMLE;

%% plot error
errCubMLE = abs(exactInteg - muhat);
plotCubatureError(dim, nvec, errCubMLE, out.ErrBd, fName, BernPolyOrder, out.ptransform, fullPath)

fprintf('done');


end
