%% Pricing Options with Control Variates
% This is a test

%% Initialization
% First we set up the basic common praramters for our examples.
clc;clearvars;
iter=1;[t,t1,t2] =deal(0); [n,n1,n2] =deal(0);

%gail.InitializeWorkspaceDisplay %initialize the workspace and the display parameters
inp.timeDim.timeVector = 1/52:1/52:52*1/52; %weekly monitoring for 1 year 
inp.assetParam.initPrice = 36; %initial stock price
inp.assetParam.interest = 0.06; %risk-free interest rate
inp.assetParam.volatility = 0.5; %volatility
inp.payoffParam.strike = 40; %strike price
inp.priceParam.absTol = 1e-2; %absolute tolerance of a nickel
inp.priceParam.relTol = 0; %zero relative tolerance
inp.priceParam.cubMethod = 'Sobol'; %Sobol sampling
inp.bmParam.assembleType = 'diff';
EuroCall = optPrice(inp); %construct an optPrice object 
opt = optPayoff(EuroCall); %make a copy
opt1 = optPayoff(EuroCall); %make a copy
opt2 = optPayoff(EuroCall); %make a copy

% american put option with european option and stockprice as c.v.s 
	opt.payoffParam = struct( ...
	'optType',{{'american'}},...
	'putCallType',{{'put'}}); 
    % single CV (euro)
	opt1.payoffParam = struct( ...
	'optType',{{'american','euro'}},...
	'putCallType',{{'put','put'}}); 
    % double CV (euro, S(T))
	opt2.payoffParam = struct( ...
	'optType',{{'american','euro','stockprice'}},...
	'putCallType',{{'put','put',''}});
% opt2.payoffParam = struct( ...
%     'strike',{{inp.payoffParam.strike, inp.payoffParam.strike, inp.payoffParam.strike}};

% make param shorter 
abstol = inp.priceParam.absTol;
reltol = inp.priceParam.relTol;
d = opt.timeDim.nSteps; 
% get exact price of the option, 0 if not known
exactf = opt.exactPrice; exactf(isnan(exactf))=0;

f =@(x) genOptPayoffs(opt,x);
f1.func =@(x) genOptPayoffs(opt1,x);
f1.cv = opt1.exactPrice; f1.cv = f1.cv(2:end); 
f2.func =@(x) genOptPayoffs(opt2,x);
f2.cv = opt2.exactPrice;f2.cv=f2.cv(2:end); 
fprintf(' abstol=%s \n ', abstol);
% begin testing no cv
for i=1:iter
	[q,out] = cubSobol_american_g(f,[zeros(1,d) ; ones(1,d)],abstol,0);
	t=t+out.time;
       	n=n+out.n;
end
fprintf('\n Results of cubSobol_american_g: \n');
fprintf('q=%8.7f  \n',q);
fprintf('err=|q-exact|=%.7f  \n',abs(q-exactf));
fprintf('avg time of cubSobol_american_g: %s \n', num2str(t/iter) ); 
fprintf('avg n of cubSobol_american_g: %s \n', num2str(n/iter) );
% begin testing single cv
for i=1:iter
    [q1,out1] = cubSobol_american_g(f1,[zeros(1,d) ; ones(1,d)], abstol,0);
    t1=t1+out1.time;
    n1=n1+out1.n;
end
fprintf('\n Results of cubSobol_american_g(single cv): \n');
fprintf('q1=%8.7f  \n',q1);
fprintf('err=|q1-exact|=%.7f  \n',abs(q1-exactf));
fprintf('avg time of cubSobol_american_g: %s \n', num2str(t1/iter) ); 
fprintf('avg n of cubSobol_american_g: %s \n', num2str(n1/iter) );
fprintf('beta %4d \n', out1.beta );

% begin testing double cv
for i=1:iter
    [q2,out2] = cubSobol_american_g(f2,[zeros(1,d) ; ones(1,d)],abstol,0);
    t2=t2+out2.time;
    n2=n2+out2.n;
end
fprintf('\n Results of cubSobol_american_g(double cv): \n');
fprintf('q2=%8.7f  \n',q2);
fprintf('err=|q2-exact|=%.7f  \n',abs(q2-exactf));
fprintf('avg time of cubSobol_american_g: %s \n', num2str(t2/iter) ); 
fprintf('avg n of cubSobol_american_g: %s \n', num2str(n2/iter) );
fprintf('beta %4d \n', out2.beta );


%
% %%
% _Author: Da Li 
