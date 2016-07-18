function [K,kvec,k0] = kernelFun(x,whKer,shape,domain)
[nx,d] = size(x);
if nargin < 4
   domain = [zeros(1,d); ones(1,d)];
   if nargin < 3;
      shape = 1;
      if nargin < 2
         whKer = 'sqExp';
      end
   end
end
K = ones(nx);
if strcmp(whKer,'sqExp')
   kvec = ones(nx,1)*(sqrt(pi)/(2*shape))^d;
   for k = 1:d;
      K = K.*exp(-(shape*bsxfun(@minus,x(:,k),x(:,k)')).^2);
      kvec = kvec.*(erf(shape*x(:,k)) + erf(shape*(1 - x(:,k))));
   end
elseif strcmp(whKer,'Mat1')
   kvec = ones(nx,1)*((2/shape)^d);
   for k = 1:d;
      tempa = shape*abs(bsxfun(@minus,x(:,k),x(:,k)'));
      K = K.*exp(-tempa).*(1 + tempa);
      tempb = shape*x(:,k);
      tempc = shape - tempb;
      kvec = kvec.*(2 - exp(-tempc).*(1+tempc/2) ...
          - exp(-tempb).*(1+tempb/2));
   end
end
    
end
