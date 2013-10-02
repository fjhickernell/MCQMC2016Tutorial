% Yizhi Zhang yzhang97@hawk.iit.edu
%
% Windows 7 Home Premium, Matlab 8.1 (R2013a)

n = 20;                                  
A = randn(n,n);
A = triu(A,1)-triu(A,1)'; % to create a skew symmetric matrix.
b = randn(n,1);
tol = 1e-8;   
maxit = 100;


tic
x = ssminresqlp(A,b,tol,maxit);
time=toc
memory

% % output:
% Enter ssminresqlp.  Min-length solution of symmetric (A-sI)x = b or min ||(A-sI)x - b||
% n      =     20   ||b||    = 4.817e+00   shift    = 0.000e+00   rtol     = 1.000e-08
% maxit  =    100   maxxnorm = 1.000e+07   Acondlim = 1.000e+15   TranCond = 1.000e+07
% precon =      0
% 
%     iter     rnorm     Arnorm   Compatible     LS        Anorm      Acond      xnorm
% P      0   4.82e+00   2.72e+01   1.00e+00   1.00e+00   0.00e+00   1.00e+00   0.00e+00
%        1   4.82e+00   2.72e+01   1.00e+00   4.26e-01   5.64e+00   1.00e+00   0.00e+00
%        2   4.36e+00   3.47e+01   6.34e-01   3.71e-01   1.32e+01   2.35e+00   1.55e-01
%        3   4.36e+00   3.47e+01   5.35e-01   2.98e-01   2.15e+01   3.81e+00   1.55e-01
%        4   4.16e+00   4.02e+01   3.87e-01   2.75e-01   2.67e+01   4.74e+00   2.22e-01
%        5   4.16e+00   4.02e+01   3.30e-01   2.38e-01   3.51e+01   6.22e+00   2.22e-01
%        6   4.04e+00   4.47e+01   2.59e-01   2.26e-01   4.06e+01   7.20e+00   2.66e-01
%        7   4.04e+00   4.47e+01   2.27e-01   2.02e-01   4.90e+01   8.68e+00   2.66e-01
%        8   3.96e+00   4.85e+01   1.88e-01   1.95e-01   5.46e+01   9.68e+00   2.97e-01
%        9   3.96e+00   4.85e+01   1.68e-01   1.78e-01   6.30e+01   1.12e+01   2.97e-01
%       10   3.89e+00   5.19e+01   1.45e-01   1.73e-01   6.88e+01   1.22e+01   3.22e-01
%       11   3.89e+00   5.19e+01   1.31e-01   1.61e-01   7.72e+01   1.37e+01   3.22e-01
%       12   3.84e+00   5.50e+01   1.16e-01   1.57e-01   8.30e+01   1.47e+01   3.42e-01
%       13   3.84e+00   5.50e+01   1.07e-01   1.47e-01   9.14e+01   1.62e+01   3.42e-01
%       14   3.80e+00   5.78e+01   9.58e-02   1.44e-01   9.73e+01   1.72e+01   3.59e-01
%       15   3.80e+00   5.78e+01   8.90e-02   1.36e-01   1.06e+02   1.87e+01   3.59e-01
%       16   3.77e+00   6.04e+01   8.11e-02   1.34e-01   1.12e+02   1.98e+01   3.73e-01
%       17   3.77e+00   6.04e+01   7.60e-02   1.27e-01   1.20e+02   2.13e+01   3.73e-01
%       18   3.74e+00   6.29e+01   6.99e-02   1.25e-01   1.26e+02   2.23e+01   3.86e-01
%       19   3.74e+00   6.29e+01   6.59e-02   1.20e-01   1.34e+02   2.38e+01   3.86e-01
%       20   3.71e+00   6.52e+01   6.12e-02   1.18e-01   1.40e+02   2.49e+01   3.97e-01
% 
%     iter     rnorm     Arnorm   Compatible     LS        Anorm      Acond      xnorm
%       21   3.71e+00   6.52e+01   5.80e-02   1.14e-01   1.49e+02   2.64e+01   3.97e-01
%       22   3.68e+00   6.74e+01   5.43e-02   1.12e-01   1.55e+02   2.74e+01   4.07e-01
%       23   3.68e+00   6.74e+01   5.17e-02   1.08e-01   1.63e+02   2.89e+01   4.07e-01
%       24   3.66e+00   6.94e+01   4.86e-02   1.07e-01   1.69e+02   3.00e+01   4.17e-01
%       25   3.66e+00   6.94e+01   4.65e-02   1.03e-01   1.78e+02   3.15e+01   4.17e-01
%       26   3.64e+00   7.14e+01   4.39e-02   1.02e-01   1.84e+02   3.26e+01   4.25e-01
%       27   3.64e+00   7.14e+01   4.21e-02   9.89e-02   1.92e+02   3.40e+01   4.25e-01
%       28   3.63e+00   7.33e+01   4.00e-02   9.79e-02   1.98e+02   3.51e+01   4.33e-01
%       29   3.63e+00   7.33e+01   3.85e-02   9.51e-02   2.07e+02   3.66e+01   4.33e-01
%       30   3.61e+00   7.51e+01   3.67e-02   9.42e-02   2.13e+02   3.77e+01   4.40e-01
%       31   3.61e+00   7.51e+01   3.53e-02   9.16e-02   2.21e+02   3.92e+01   4.40e-01
%       32   3.59e+00   7.69e+01   3.38e-02   9.08e-02   2.27e+02   4.03e+01   4.47e-01
%       33   3.59e+00   7.69e+01   3.26e-02   8.85e-02   2.36e+02   4.18e+01   4.47e-01
%       34   3.58e+00   7.86e+01   3.13e-02   8.78e-02   2.42e+02   4.29e+01   4.53e-01
%       35   3.58e+00   7.86e+01   3.03e-02   8.57e-02   2.50e+02   4.43e+01   4.53e-01
%       36   3.57e+00   8.03e+01   2.91e-02   8.50e-02   2.56e+02   4.54e+01   4.59e-01
%       37   3.57e+00   8.03e+01   2.82e-02   8.31e-02   2.65e+02   4.69e+01   4.59e-01
%       38   3.55e+00   8.18e+01   2.72e-02   8.25e-02   2.71e+02   4.80e+01   4.65e-01
%       39   3.55e+00   8.18e+01   2.64e-02   8.07e-02   2.79e+02   4.95e+01   4.65e-01
%       40   3.54e+00   8.34e+01   2.55e-02   8.02e-02   2.85e+02   5.06e+01   4.70e-01
% 
%     iter     rnorm     Arnorm   Compatible     LS        Anorm      Acond      xnorm
%       41   3.54e+00   8.34e+01   2.48e-02   7.85e-02   2.94e+02   5.21e+01   4.70e-01
%       42   3.53e+00   8.49e+01   2.40e-02   7.80e-02   3.00e+02   5.32e+01   4.75e-01
%       43   3.53e+00   8.49e+01   2.34e-02   7.64e-02   3.08e+02   5.47e+01   4.75e-01
%       44   3.52e+00   8.64e+01   2.26e-02   7.60e-02   3.15e+02   5.58e+01   4.80e-01
%       45   3.52e+00   8.64e+01   2.21e-02   7.45e-02   3.23e+02   5.72e+01   4.80e-01
%       46   3.51e+00   8.78e+01   2.14e-02   7.41e-02   3.29e+02   5.83e+01   4.84e-01
%       47   3.51e+00   8.78e+01   2.09e-02   7.27e-02   3.37e+02   5.98e+01   4.84e-01
%       48   3.50e+00   8.92e+01   2.03e-02   7.23e-02   3.44e+02   6.09e+01   4.89e-01
%       49   3.50e+00   8.92e+01   1.98e-02   7.11e-02   3.52e+02   6.24e+01   4.89e-01
%       50   3.49e+00   9.06e+01   1.93e-02   7.07e-02   3.58e+02   6.35e+01   4.93e-01
%       51   3.49e+00   9.06e+01   1.88e-02   6.95e-02   3.67e+02   6.50e+01   4.93e-01
%       52   3.48e+00   9.19e+01   1.83e-02   6.92e-02   3.73e+02   6.61e+01   4.97e-01
%       53   3.48e+00   9.19e+01   1.80e-02   6.80e-02   3.81e+02   6.76e+01   4.97e-01
%       54   3.48e+00   9.32e+01   1.75e-02   6.77e-02   3.88e+02   6.87e+01   5.00e-01
%       55   3.48e+00   9.32e+01   1.71e-02   6.67e-02   3.96e+02   7.02e+01   5.00e-01
%       56   3.47e+00   9.45e+01   1.67e-02   6.63e-02   4.02e+02   7.13e+01   5.04e-01
%       57   3.47e+00   9.45e+01   1.64e-02   6.53e-02   4.10e+02   7.28e+01   5.04e-01
%       58   3.46e+00   9.57e+01   1.60e-02   6.51e-02   4.17e+02   7.39e+01   5.07e-01
%       59   3.46e+00   9.57e+01   1.57e-02   6.41e-02   4.25e+02   7.54e+01   5.07e-01
%       60   3.45e+00   9.70e+01   1.53e-02   6.38e-02   4.31e+02   7.65e+01   5.11e-01
% 
%     iter     rnorm     Arnorm   Compatible     LS        Anorm      Acond      xnorm
%       61   3.45e+00   9.70e+01   1.51e-02   6.29e-02   4.40e+02   7.80e+01   5.11e-01
%       62   3.45e+00   9.82e+01   1.47e-02   6.27e-02   4.46e+02   7.91e+01   5.14e-01
%       63   3.45e+00   9.82e+01   1.45e-02   6.18e-02   4.54e+02   8.05e+01   5.14e-01
%       64   3.44e+00   9.93e+01   1.42e-02   6.16e-02   4.61e+02   8.17e+01   5.17e-01
%       65   3.44e+00   9.93e+01   1.39e-02   6.07e-02   4.69e+02   8.31e+01   5.17e-01
%       66   3.43e+00   1.01e+02   1.36e-02   6.05e-02   4.75e+02   8.43e+01   5.20e-01
%       67   3.43e+00   1.01e+02   1.34e-02   5.97e-02   4.84e+02   8.57e+01   5.20e-01
%       68   3.43e+00   1.02e+02   1.31e-02   5.95e-02   4.90e+02   8.69e+01   5.23e-01
%       69   3.43e+00   1.02e+02   1.29e-02   5.87e-02   4.98e+02   8.83e+01   5.23e-01
%       70   3.42e+00   1.03e+02   1.27e-02   5.85e-02   5.05e+02   8.95e+01   5.26e-01
%       71   3.42e+00   1.03e+02   1.25e-02   5.78e-02   5.13e+02   9.09e+01   5.26e-01
%       72   3.42e+00   1.04e+02   1.22e-02   5.76e-02   5.19e+02   9.21e+01   5.29e-01
%       73   3.42e+00   1.04e+02   1.20e-02   5.69e-02   5.28e+02   9.35e+01   5.29e-01
%       74   3.41e+00   1.05e+02   1.18e-02   5.67e-02   5.34e+02   9.47e+01   5.31e-01
%       75   3.41e+00   1.05e+02   1.16e-02   5.61e-02   5.42e+02   9.61e+01   5.31e-01
%       76   3.41e+00   1.06e+02   1.14e-02   5.59e-02   5.49e+02   9.73e+01   5.34e-01
%       77   3.41e+00   1.06e+02   1.13e-02   5.53e-02   5.57e+02   9.87e+01   5.34e-01
%       78   3.40e+00   1.07e+02   1.11e-02   5.51e-02   5.63e+02   9.99e+01   5.37e-01
%       79   3.40e+00   1.07e+02   1.09e-02   5.45e-02   5.72e+02   1.01e+02   5.37e-01
%       80   3.40e+00   1.08e+02   1.07e-02   5.43e-02   5.78e+02   1.02e+02   5.39e-01
% 
%     iter     rnorm     Arnorm   Compatible     LS        Anorm      Acond      xnorm
%       81   3.40e+00   1.08e+02   1.06e-02   5.37e-02   5.86e+02   1.04e+02   5.39e-01
%       82   3.39e+00   1.09e+02   1.04e-02   5.36e-02   5.93e+02   1.05e+02   5.41e-01
%       83   3.39e+00   1.09e+02   1.03e-02   5.30e-02   6.01e+02   1.07e+02   5.41e-01
%       84   3.39e+00   1.10e+02   1.01e-02   5.28e-02   6.07e+02   1.08e+02   5.44e-01
%       85   3.39e+00   1.10e+02   9.97e-03   5.23e-02   6.16e+02   1.09e+02   5.44e-01
%       86   3.38e+00   1.11e+02   9.82e-03   5.22e-02   6.22e+02   1.10e+02   5.46e-01
%       87   3.38e+00   1.11e+02   9.69e-03   5.16e-02   6.30e+02   1.12e+02   5.46e-01
%       88   3.38e+00   1.12e+02   9.55e-03   5.15e-02   6.37e+02   1.13e+02   5.48e-01
%       89   3.38e+00   1.12e+02   9.42e-03   5.10e-02   6.45e+02   1.14e+02   5.48e-01
%       90   3.37e+00   1.13e+02   9.28e-03   5.08e-02   6.51e+02   1.15e+02   5.50e-01
%       91   3.37e+00   1.13e+02   9.17e-03   5.04e-02   6.60e+02   1.17e+02   5.50e-01
%       92   3.37e+00   1.14e+02   9.04e-03   5.02e-02   6.66e+02   1.18e+02   5.52e-01
%       93   3.37e+00   1.14e+02   8.93e-03   4.97e-02   6.74e+02   1.20e+02   5.52e-01
%       94   3.36e+00   1.15e+02   8.80e-03   4.96e-02   6.81e+02   1.21e+02   5.54e-01
%       95   3.36e+00   1.15e+02   8.70e-03   4.92e-02   6.89e+02   1.22e+02   5.54e-01
%       96   3.36e+00   1.16e+02   8.58e-03   4.90e-02   6.95e+02   1.23e+02   5.56e-01
%       97   3.36e+00   1.16e+02   8.48e-03   4.86e-02   7.04e+02   1.25e+02   5.56e-01
%       98   3.36e+00   1.17e+02   8.36e-03   4.85e-02   7.10e+02   1.26e+02   5.58e-01
%       99   3.36e+00   1.17e+02   8.27e-03   4.80e-02   7.18e+02   1.27e+02   5.58e-01
%      100   6.64e+00   4.19e+01D  2.95e-02   8.71e-03D  7.25e+02   1.28e+02   3.04e-01
% 
% 
% Exit ssminresqlp.   flag  =      8   The iteration limit was reached                        
% Exit ssminresqlp.   iter  =    100   (SS-MINRES      0, SS-MINRES-QLP    100)
% Exit ssminresqlp.   rnorm =  6.6420e+00     rnorm  direct =  6.6420e+00
% Exit ssminresqlp.                           Arnorm direct =  4.1933e+01
% Exit ssminresqlp.   xnorm =  3.0429e-01     xnorm  direct =  3.0429e-01
% Exit ssminresqlp.   Anorm =  7.2481e+02     Acond         =  1.2849e+02
% time =
%     0.0484
% Maximum possible array:      4450 MB (4.666e+09 bytes) *
% Memory available for all arrays:      4450 MB (4.666e+09 bytes) *
% Memory used by MATLAB:       900 MB (9.435e+08 bytes)
% Physical Memory (RAM):      4055 MB (4.252e+09 bytes)
% 
% *  Limited by System Memory (physical + swap file) available.


% Suggestion for improvement:  
% 1.  In file shminresqlp.m, the function name is incorrect. It appears
%"minresqlp(A,b,rtol,maxit,M,shift,maxxnorm,Acondlim,TranCond,show)".
%However it should be
%"shminresqlp(A,b,rtol,maxit,M,shift,maxxnorm,Acondlim,TranCond,show)". 
%
% 2. Moreover, the first line in of the documentation (description) is not
% correct.
% 
% 3. For the ssminresqlp.m, if we specify the input "condition number",
% line 256 is operated, which leads to failure in function minresxxxM, line
% 665. I have tried with the examples in side. They failed, too.