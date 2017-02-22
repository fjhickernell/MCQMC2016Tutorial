function [val] = MLEKernel(shape,x,y,whKer,domain,BernPolynX,BernPolynOrder)
% The function MLEKERNEL is the loss function to be minimized to obtain the
% shape parameter |shape| that defines the kernel.  It uses data |(x,y)|
% and the function |kernelFun|.
if nargin < 4
    whKer = 'Mat1';
end
nx = size(x,1);

K = kernelFun(x,whKer,shape,domain,BernPolynX,BernPolynOrder);
if strcmp(whKer,'Fourier')
    cn = K; % only first row is given
    
    eigval = abs(fft((cn-1)'))'/nx;   % this is done to improve accuracy, to reduce zero values
    eigval(1) = 1;
    if any((eigval==0))
        fprintf('zero_eig %d, ', sum(eigval==0));
    end
    
    val1 = sum(log( eigval(eigval~=0) ))/nx ;
    
    ffy = abs(fft(y-1))'/(nx);   % this is done to improve accuracy, to reduce zero values
    ffy(1) = 1;
    temp = ((ffy(eigval~=0).^2)./abs(eigval(eigval~=0)));
    val2 = log(sum(temp));
    
    val = val1 + val2;
    if isnan(val)
        fprintf('nan ');
    end
    val = real(val);
    
else
    [eigvec,eigval] = eig(K,'vector');
    Vty = eigvec'*y;
    val = sum(log(eigval))/nx + log(Vty'*(Vty./eigval));
    %val = sum(log(eigval))/nx + Vty'*(Vty./eigval); %original code
end

end


