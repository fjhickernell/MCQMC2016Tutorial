function val = MLEKernel(shape,x,y,whKer,domain,BernPolynX,BernPolynOrder)
% The function MLEKERNEL is the loss function to be minimized to obtain the
% shape parameter |shape| that defines the kernel.  It uses data |(x,y)|
% and the function |kernelFun|.
if nargin < 4
   whKer = 'Mat1';
end
nx = size(x,1);
n1 = nx-1;
K = kernelFun(x,whKer,shape,domain,BernPolynX,BernPolynOrder);
if strcmp(whKer,'Fourier')
    cn = K(1,:);
    wj = exp(2*pi*1i*(0:n1)'/nx);
    V = fliplr(vander(wj));
    eigval = V*cn';
    eigvec = (1/sqrt(nx))*V;
else
    [eigvec,eigval] = eig(K,'vector');
end
Vty = eigvec'*y;
val = sum(log(eigval))/nx + Vty'*(Vty./eigval);
val = real(val);
end
