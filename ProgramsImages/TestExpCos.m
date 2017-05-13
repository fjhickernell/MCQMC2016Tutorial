%% Test out MLE
function [nvec,muhat,aMLE,errCubMLE,out] = TestExpCos(dim,BernPolyOrder,ptransform,figSavePath)

nvec = 2.^(8:20);

domain = repmat([0;1],[1,dim]);

f = @(x) exp(sum(cos(2*pi*x), 2));
exactInteg = besseli(0,1)^dim;
fName = 'Exp(cos)';
%figSavePath = '/home/jagadees/MyWriteup/Apr1stweek/';

whSample = 'Lattice1';
whKer = 'Fourier';
powerFuncMethod = 'Cauchy';

fullPath = strcat(figSavePath,'/',fName,'/',ptransform,'/');
if exist(fullPath,'dir')==false
    mkdir(fullPath);
end

if 0
    [muhat,out] = cubMLE(f,nvec,domain, whSample,whKer,powerFuncMethod,BernPolyOrder,...
        ptransform,fName,fullPath);
    aMLE = out.aMLE;
    ErrBd = out.ErrBd;
else
    absTol = 1E-15;
    relTol = 0;
    order = BernPolyOrder;
    testAll = true;
    [muhatFinal,out]=cubMLELattice(f,dim,absTol,relTol,order,ptransform,testAll,figSavePath,fName);
    nvec = 2.^out.mvec;
    muhat = out.muhatAll;
    ErrBd = out.ErrBdAll;
end

%% plot error
errCubMLE = abs(exactInteg - muhat);
plotCubatureError(dim, nvec, errCubMLE, ErrBd, fName, BernPolyOrder, ptransform, fullPath)

fprintf('done');


end
