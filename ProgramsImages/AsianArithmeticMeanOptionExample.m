%% Generate Examples of Asian Arithmetic Mean Option Pricing

%gail.InitializeWorkspaceDisplay %clean up 
clearvars
format long

if exist('AsianCallExampleAllData.mat','file')
   load AsianCallExampleAllData
   ArchEuroCall = AsianCall;
   ArchAsianCall = AsianCall;
   ArchAsianCallCV = AsianCall;
   Archnvec = nvec;
   ArchabsTolGold = absTolGold;
   ArchnGoldRep = nGoldRep;
   ArchAbsTol = absTol;
   ArchRelTol = relTol;
   ArchnRepAuto = nRepAuto;
   ArchnRep = nRep;
end
nvec = 2.^(7:17)';
nmax = max(nvec);
nRep = 100;
nlarge = nmax*2;
nn = numel(nvec);
alpha = 0.1;


%% Parameters for the Asian option
absTol = 1e-4;
relTol = 0;
inp.timeDim.timeVector = 1/52:1/52:1/4; %weekly monitoring for one quarter 
%inp.timeDim.timeVector = 1/4:1/4:1/4; %weekly monitoring for one quarter 
%inp.timeDim.timeVector = 1/8:1/8:1/4; %weekly monitoring for one quarter 
inp.assetParam.initPrice = 100; %initial stock price
inp.assetParam.interest = 0.01; %risk-free interest rate
inp.assetParam.volatility = 0.5; %volatility
inp.payoffParam.strike = 100; %strike price
inp.priceParam.absTol = absTol; %absolute tolerance of a penny
inp.priceParam.relTol = relTol; %zero relative tolerance
inp.priceParam.cubMethod = 'Sobol'; %Sobol sampling
inp.bmParam.assembleType = 'PCA';
inp.payoffParam.putCallType = {'call'};

%% Construct some different options
EuroCall = optPrice(inp); %construct a European optPrice object
AsianCall = optPrice(EuroCall); %construct an Asian optPrice object
AsianCall.payoffParam = struct( ...
	'optType',{{'amean'}},...
	'putCallType',{{'call'}});
AsianCallCV = optPrice(EuroCall); %construct an Asian and European optPayoff object for CV
AsianCallCV.payoffParam = struct( ...
	'optType',{{'amean','gmean'}},...
	'putCallType',{{'call','call'}});

%% Construct a very accurate answer
compGold = true;
absTolGold = 1e-6;
nGoldRep = 100;
if exist('callPriceExact','var') && ...
   absTolGold == ArchabsTolGold && nGoldRep == ArchnGoldRep && ...
   all(ArchAsianCall.timeDim.timeVector == AsianCall.timeDim.timeVector) && ...
   ArchAsianCall.assetParam.initPrice == AsianCall.assetParam.initPrice && ... %initial stock price
   ArchAsianCall.payoffParam.strike == ArchAsianCall.payoffParam.strike, %strike price   
   compGold = false;
   disp('Already have gold standard Asian Call')
end
 
if compGold
   fCV.func = @(x) genOptPayoffs(AsianCallCV,x);
   fCV.cv = AsianCallCV.exactPrice(2:end); 
   d = AsianCallCV.timeDim.nSteps;
   callPriceGold(nGoldRep,1) = 0;  
   tic
   for ii = 1:nGoldRep
      gail.TakeNote(ii,1) %print out every 10th ii
      callPriceGold(ii) = ...
         cubSobol_g(fCV,[zeros(1,d); ones(1,d)],'uniform',absTolGold,relTol);
%       x = net(scramble(sobolset(AsianCall.timeDim.nSteps), ...
%          'MatousekAffineOwen'),nGold);
%       payoffAsianEuro = genOptPayoffs(AsianCall,x);
%       callPriceGold(ii) = mean(payoffAsianEuro(:,1));
%       payoffAsianEuro = genOptPayoffs(AsianCallCV,x);
%       temp = payoffAsianEuro(:,1) - AsianCallCV.exactPrice(2) + payoffAsianEuro(:,2);
%       callPriceGold(ii) = mean(temp);
   end
   callPriceExact = mean(callPriceGold);
   toc
end
disp(['mu  = ' num2str(callPriceExact,15) ' +/- ' num2str(2*std(callPriceGold),10)])

%% IID sampling
AsianCallIID = optPayoff(AsianCall);
AsianCallIID.wnParam = struct('sampleKind','IID','xDistrib','Gaussian');
AsianCallIID.bmParam.assembleType = 'diff';
AsianCallIID.inputType = 'n';
compIID = true;
if exist('AsianCallExampleAllData.mat','file')
   if exist('muAsianCallIID','var') && numel(nvec) == numel(Archnvec) && ...
         nRep == ArchnRep
      if all(nvec == Archnvec)
         compIID = false;
         disp('Already have IID Asian Call')
      end
   end
end
if compIID
   tic
   muAsianCallIID = zeros(nn,nRep);
   for i = 1:nRep
      temp = cumsum(genOptPayoffs(AsianCallIID,nmax));
      muAsianCallIID(:,i) = temp(nvec)./nvec; 
   end
   errvecAsianCallIID = abs(callPriceExact - muAsianCallIID);
   errmedAsianCallIID = median(errvecAsianCallIID,2);
   errtopAsianCallIID = quantile(errvecAsianCallIID,1-alpha,2);
   toc
end

%% Unscrambled Sobol sampling
compUSobol = true;
if exist('AsianCallExampleAllData.mat','file')
   if exist('muAsianCallUSobol','var') && numel(nvec) == numel(Archnvec) && ...
         nRep == ArchnRep
      if all(nvec == Archnvec)
         compUSobol = false;
         disp('Already have unscrambled Sobol Asian Call')
      end
   end
end
if compUSobol
   tic 
   muAsianCallUSobol = zeros(nn,1);
     x = net(sobolset(AsianCall.timeDim.nSteps),nmax);
     temp = cumsum(genOptPayoffs(AsianCall,x));
     muAsianCallUSobol(:) = temp(nvec)./nvec; 
   errvecAsianCallUSobol = abs(callPriceExact - muAsianCallUSobol);
   toc
end

%% Scrambled Sobol sampling
compSobol = true;
if exist('AsianCallExampleAllData.mat','file')
   if exist('muAsianCallSobol','var') && numel(nvec) == numel(Archnvec) && ...
         nRep == ArchnRep
      if all(nvec == Archnvec)
         compSobol = false;
         disp('Already have scrambled Sobol Asian Call')
      end
   end
end
if compSobol
   tic 
   muAsianCallSobol = zeros(nn,nRep);
   for i = 1:nRep
      x = net(scramble(sobolset(AsianCall.timeDim.nSteps), ...
         'MatousekAffineOwen'),nmax);
      temp = cumsum(genOptPayoffs(AsianCall,x));
      muAsianCallSobol(:,i) = temp(nvec)./nvec; 
   end
   errvecAsianCallSobol = abs(callPriceExact - muAsianCallSobol);
   errmedAsianCallSobol = median(errvecAsianCallSobol,2);
   errtopAsianCallSobol = quantile(errvecAsianCallSobol,1-alpha,2);
   toc
end

%% Scrambled Sobol sampling with time differencing
AsianCallDiff = optPayoff(AsianCall);
AsianCallDiff.bmParam.assembleType = 'diff';
compSobolDiff = true;
if exist('AsianCallExampleAllData.mat','file')
   if exist('muAsianCallDiff','var') && numel(nvec) == numel(Archnvec) && ...
         nRep == ArchnRep
      if all(nvec == Archnvec)
         compSobolDiff = false;
         disp('Already have scrambled Sobol Asian Call w/ time differencing')
      end
   end
end
if compSobolDiff
   tic 
   muAsianCallDiff = zeros(nn,nRep);
   for i = 1:nRep
      x = net(scramble(sobolset(AsianCall.timeDim.nSteps), ...
         'MatousekAffineOwen'),nmax);
      temp = cumsum(genOptPayoffs(AsianCallDiff,x));
      muAsianCallDiff(:,i) = temp(nvec)./nvec; 
   end
   errvecAsianCallDiffSobol = abs(callPriceExact - muAsianCallDiff);
   errmedAsianCallDiffSobol = median(errvecAsianCallDiffSobol,2);
   errtopAsianCallDiffSobol = quantile(errvecAsianCallDiffSobol,1-alpha,2);
   toc
end

%% Try automatic cubature
nRepAuto = 100;
compCallAuto = true;
if exist('AsianCallExampleAllData.mat','file')
   if exist('muAsianCallAuto','var') && ...
      numel(nRepAuto) == numel(ArchnRepAuto) && ...
      ArchAbsTol == absTol && ...
      ArchRelTol == relTol 
      compCallAuto = false;
      disp('Already have automatic scrambled Sobol Asian Call')
   end
end
if compCallAuto
   f = @(x) genOptPayoffs(AsianCall,x);
   d = AsianCall.timeDim.nSteps;
   muAsianCallAuto(nRepAuto,1) = 0;
   tic
   for i =  1:nRepAuto
      gail.TakeNote(i,10)
      [muAsianCallAuto(i),outCall(i)] = ...
         cubSobol_g(f,[zeros(1,d); ones(1,d)],'uniform',absTol,relTol);
   end
   errvecAsianCallAuto = abs(callPriceExact - muAsianCallAuto);
   errmedAsianCallAuto = median(errvecAsianCallAuto,2);
   errtopAsianCallAuto = quantile(errvecAsianCallAuto,1-alpha,2);
   rangeAsianCallAuto = range(muAsianCallAuto);
   toc
end

%% Try automatic control variates
compCallCVAuto = true;
if exist('AsianCallExampleAllData.mat','file')
   if exist('muAsianCallCVAuto','var') && ...
      numel(nRepAuto) == numel(ArchnRepAuto) && ...
      ArchAbsTol == absTol && ...
      ArchRelTol == relTol 
      compCallCVAuto = false;
      disp('Already have automatic scrambled Sobol Asian Call with control variates')
   end
end
if compCallCVAuto
   fCV.func = @(x) genOptPayoffs(AsianCallCV,x);
   fCV.cv = AsianCallCV.exactPrice(2:end); 
   d = AsianCallCV.timeDim.nSteps;
   muAsianCallCVAuto(nRepAuto,1) = 0;
   tic
   for i =  1:nRepAuto
      gail.TakeNote(i,10)
      [muAsianCallCVAuto(i),outCallCV(i)] = ...
         cubSobol_g(fCV,[zeros(1,d); ones(1,d)],'uniform',absTol,relTol);
   end
   errvecAsianCallCVAuto = abs(callPriceExact - muAsianCallCVAuto);
   errmedAsianCallCVAutoo = median(errvecAsianCallCVAuto,2);
   errtopAsianCallCVAuto = quantile(errvecAsianCallCVAuto,1-alpha,2);
   rangeAsianCallCVAuto = range(muAsianCallCVAuto);
   toc
end

%% Save output
save AsianCallExampleAllData
return

