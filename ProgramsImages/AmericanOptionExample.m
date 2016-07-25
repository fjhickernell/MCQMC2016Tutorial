%% Generate Examples of American Option Pricing

InitializeWorkspaceDisplay %clean up 
format long

nvec = 2.^(7:20)';
nlarge = nvec(end)*2;
nn = numel(nvec);
alpha = 0.1;

%% Parameters for the American option
absTol = 1e-3;
relTol = 0;
inp.timeDim.timeVector = 1/52:1/52:1/4; %weekly monitoring for one quarter 
%inp.timeDim.timeVector = 1/4:1/4:1/4; %weekly monitoring for one quarter 
%inp.timeDim.timeVector = 1/8:1/8:1/4; %weekly monitoring for one quarter 
inp.assetParam.initPrice = 100; %initial stock price
inp.assetParam.interest = 0.05; %risk-free interest rate
inp.assetParam.volatility = 0.5; %volatility
inp.payoffParam.strike = 100; %strike price
inp.priceParam.absTol = absTol; %absolute tolerance of a penny
inp.priceParam.relTol = relTol; %zero relative tolerance
inp.priceParam.cubMethod = 'Sobol'; %Sobol sampling
inp.bmParam.assembleType = 'PCA';
inp.payoffParam.putCallType = {'put'};

%% Construct some different options
EuroPut = optPrice(inp); %construct a European optPrice object
AmerPut = optPrice(EuroPut); %construct an American optPrice object
AmerPut.payoffParam = struct( ...
	'optType',{{'american'}},...
	'putCallType',{{'put'}});
nRep = 10;
AmerPrice(nRep,1) = 0;
%outAmerPut = struct([]);
%outAmerPut(nRep,1) = 0;
f = @(x) genOptPayoffs(AmerPut,x);
AmerOptFourierCoeffDecay(f,AmerPut.timeDim.nSteps)
return
for ii = 1:nRep
   [AmerPrice(ii),outAmerPut(ii)] = genOptPrice(AmerPut);
end

AmerPutCV = optPayoff(EuroPut); %construct an American and European optPayoff object for CV
AmerPutCV.payoffParam = struct( ...
	'optType',{{'american','euro'}},...
	'putCallType',{{'put','put'}});
f.func =@(x) genOptPayoffs(AmerPutCV,x);
f.cv = AmerPutCV.exactPrice(2:end); 
d = AmerPutCV.timeDim.nSteps;
n = 2^20;
y = f.func(net(scramble(sobolset(d),'MatousekAffineOwen'),n));
% figure
% plot(y(:,2),y(:,1),'.')
% xlabel('European payoffs')
% ylabel('American payoffs')
% corrmat = corr(y)
% covmat = cov(y);
% beta = covmat(1,2)/covmat(2,2)
% stdErrNotCorrectY = 2*std(y(:,1))/sqrt(n)
% stdErrNotCorrectYCV = 2*std(y(:,1) + beta * (AmerPutCV.exactPrice(2:end) - y(:,2)))/sqrt(n)
AmerPriceCV(nRep,1) = 0;
%outAmerPutCV(nRep,1) = 0;
for ii = 1:nRep
   [AmerPriceCV(ii),outAmerPutCV(ii)] ...
      = cubSobol_american_g(f, ...
      [zeros(1,AmerPutCV.timeDim.nSteps); ones(1,AmerPutCV.timeDim.nSteps)], ...
      absTol, relTol);
end
comparePrice = [AmerPrice AmerPriceCV]
range(comparePrice,1)
range(comparePrice(:))
compareN = [outAmerPut(:).nPaths outAmerPutCV(:).n]


