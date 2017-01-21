function val = MLEKernel(shape,x,y,whKer,domain,BernPolynX,BernPolynOrder)
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
    eigval = real(fft(cn));
    val = sum(log( eigval ))/nx + log(y'*ifft(fft(y)./eigval'));
else
    [eigvec,eigval] = eig(K,'vector');
    Vty = eigvec'*y;
    val = sum(log(eigval))/nx + log(Vty'*(Vty./eigval));
    %val = sum(log(eigval))/nx + Vty'*(Vty./eigval); %original code
end

end
