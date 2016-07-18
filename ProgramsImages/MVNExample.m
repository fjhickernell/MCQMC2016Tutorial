%% Generate Examples of Multivariate Normal Probabilities

InitializeWorkspaceDisplay %clean up 
format long

nvec = 2.^(7:20)';
nlarge = nvec(end)*2;
nn = numel(nvec);
C = [4 1 1; 0 1 0.5; 0 0 0.25];
Cov = C'*C
mu = 0;
a = [-6 -2 -2];
b = [5 2 1];
nRep = 100;
alpha = 0.1;

if exist('MVNProbExampleAllData.mat','file')
   load MVNProbExampleAllData
   MVNProbBestArch = MVNProbBest;
   nRepGoldArch = nRepGold;
   nRepArch = nRep;
   MVNProbIIDGnArch = MVNProbIIDGn;
   MVNProbSobolGnArch = MVNProbSobolGn;
   MVNProbuSobolGnArch = MVNProbuSobolGn;
   MVNProbOptSobolGnArch = MVNProbMLESobolGn;
end


%% First compute a high accuracy answer
nGold = 2^27;
nRepGold = nRep;
MVNProbBest = multivarGauss('a',a,'b',b,'Cov',Cov,'n',nGold, ...
   'errMeth','n','cubMeth','Sobol','intMeth','Genz');
compGold = true;
if exist('MVNProbExampleAllData.mat','file')
   if sameProblem(MVNProbBest,MVNProbBestArch) && ...
      nRepGoldArch == nRepGold
      disp('Already have gold standard answer')
      compGold = false;
   end
end
if compGold
   disp('(Re-)computing gold standard answer')
   muBestvec = zeros(1,nRepGold);
   tic 
   for i = 1:nRepGold
      i
      muBestvec(1,i) = compProb(MVNProbBest); 
   end
   toc
   muBest = mean(muBestvec);
end
disp(['mu  = ' num2str(muBest,15) ' +/- ' num2str(2*std(muBestvec),10)])

%% IID sampling
MVNProbIIDGn = multivarGauss('a',a,'b',b,'Cov',Cov,'n',nvec, ...
   'errMeth','n','cubMeth','IID','intMeth','Genz');
compIID = true;
if exist('MVNProbExampleAllData.mat','file')
   if sameProblem(MVNProbIIDGn,MVNProbIIDGnArch) && ...
      all(nRep == nRepArch)
      disp('Already have IID answer')
      compIID = false;
   end
end
if compIID
   tic
   muMVNProbIIDGn = zeros(nn,nRep);
   for i = 1:nRep
      muMVNProbIIDGn(:,i) = compProb(MVNProbIIDGn); 
   end
   errvecMVNProbIIDGn = abs(muBest - muMVNProbIIDGn);
   errmedMVNProbIIDGn = median(errvecMVNProbIIDGn,2);
   errtopMVNProbIIDGn = quantile(errvecMVNProbIIDGn,1-alpha,2);
   toc
end


%% Scrambled Sobol sampling
MVNProbSobolGn = multivarGauss('a',a,'b',b,'Cov',Cov,'n',nvec, ...
   'errMeth','n','cubMeth','Sobol','intMeth','Genz');
compSobol = true;
if exist('MVNProbExampleAllData.mat','file')
   if sameProblem(MVNProbSobolGn,MVNProbSobolGnArch)
      disp('Already have Scrambled Sobol answer')
      compSobol = false;
   end
end
if compSobol
   tic 
   muMVNProbSobolGn = zeros(nn,nRep);
   for i = 1:nRep
      muMVNProbSobolGn(:,i) = compProb(MVNProbSobolGn); 
   end
   errvecMVNProbSobolGn = abs(muBest - muMVNProbSobolGn);
   errmedMVNProbSobolGn = median(errvecMVNProbSobolGn,2);
   errtopMVNProbSobolGn = quantile(errvecMVNProbSobolGn,1-alpha,2);
   toc
end

%% Scrambled Sobol sampling for Affine transformation
MVNProbSobolAn = multivarGauss('a',a,'b',b,'Cov',Cov,'n',nvec, ...
   'errMeth','n','cubMeth','Sobol','intMeth','aff');
compSobolaff = true;
if exist('MVNProbExampleAllData.mat','file')
   if sameProblem(MVNProbSobolGn,MVNProbSobolGnArch)
      disp('Already have affine Scrambled Sobol answer')
      compSobolaff = false;
   end
end
if compSobolaff
   tic 
   muMVNProbSobolAn = zeros(nn,nRep);
   for i = 1:nRep
      muMVNProbSobolAn(:,i) = compProb(MVNProbSobolAn); 
   end
   errvecMVNProbSobolAn = abs(muBest - muMVNProbSobolAn);
   errmedMVNProbSobolAn = median(errvecMVNProbSobolAn,2);
   errtopMVNProbSobolAn = quantile(errvecMVNProbSobolAn,1-alpha,2);
   toc
end

%% Unscrambled Sobol sampling
MVNProbuSobolGn = multivarGauss('a',a,'b',b,'Cov',Cov,'n',nvec, ...
   'errMeth','n','cubMeth','uSobol','intMeth','Genz');
compuSobol = true;
if exist('MVNProbExampleData.mat','file')
   if sameProblem(MVNProbuSobolGn,MVNProbuSobolGnArch)
      disp('Already have unscrambled Sobol answer')
      compuSobol = false;
   end
end
if compuSobol
   tic
   muMVNProbuSobolGn = compProb(MVNProbuSobolGn); 
   errMVNProbuSobolGn = abs(muBest - muMVNProbuSobolGn);
   toc
end

%% Try MLE Bayseian cubature
nvecMLE = 2.^(7:11)';
nnMLE = numel(nvecMLE);
MVNProbMLESobolGn = multivarGauss('a',a,'b',b,'Cov',Cov,'n',nvecMLE, ...
   'errMeth','n','cubMeth','SobolOpt','intMeth','Genz');
compMLESobol = true;
% if exist('MVNProbExampleData.mat','file')
%    if sameProblem(MVNProbMLESobolGn,MVNProbMLESobolGnArch)
%       disp('Already have MLE Sobol answer')
%       compOptSobol = false;
%    end
% end
if compMLESobol
   tic
   muMVNProbMLESobolGn = zeros(nnMLE,nRep);
   tic 
   for i = 1:nRep
      if i/1 == floor(i/1), i, end
      muMVNProbMLESobolGn(:,i) = compProb(MVNProbMLESobolGn); 
   end
   errvecMVNProbMLESobolGn = abs(muBest - muMVNProbMLESobolGn);
   errmedMVNProbMLESobolGn = median(errvecMVNProbMLESobolGn,2);
   errtopMVNProbMLESobolGn = quantile(errvecMVNProbMLESobolGn,1-alpha,2);
   toc
end

%% Save output
save MVNProbExampleAllData.mat

