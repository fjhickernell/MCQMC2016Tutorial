%% Test out MLE
function [muhat,aMLE,err,out] = TestOutMLE(nvec,domain,f,exactInteg)

if nargin < 3
   f = @(x) cos(x(:,1) + x(:,2)); exactInteg = -cos(2) + 2*cos(1) - 1;
   if nargin < 2
      domain = [0 0;1 1];
      if nargin < 1
         nvec = 2.^6;
      end
   end
end
[muhat,out] = cubMLE(f,nvec,domain);
aMLE = out.aMLE;
err = abs(exactInteg - muhat);
end