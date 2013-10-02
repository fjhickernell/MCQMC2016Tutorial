%% cubMC_g
% |Monte Carlo method to evaluate a multidimentional integral to
% within a specified absolute error tolerance with guaranteed uncertainy
% within alpha.|
%
%% Syntax
% [Q,out_param] = *cubMC_g*(f)
%
% Q = *cubMC_g*(f,interval,measure,abstol,alpha,n_sigma,fudge)
%
% Q = *cubMC_g*(f,interval,'measure','uniform','abstol',abstol,'alpha',alpha,
%              'n_sigma',n_sigma,fudge',fudge)
%
% Q = *cubMC_g*(f,interval,,in_param)
%
%% Description
% [Q,out_param] = *cubMC_g*(f) |estimates the integral with integrand f to
% within the absolute error tolerance 1e-2 and with guaranteed
% uncertainty alpha within 1%. Input f a function handle. The function
% Y=f(X) should accept a vector argument X and return a vector result Y,
% the integrand evaluated at each element of X.|
%
% Q = *cubMC_g*(f,interval,measure,abstol,alpha,n_sigma,fudge) |estimates the
% integral with integrand f to within an absolute error tolerance abstol
% with guaranteed uncertainty within alpha using ordered parameter input
% interval, measure, tolerence, uncertainty, n_sigma and fudge factor.|
%
% Q = *cubMC_g*(f,interval,'measure','uniform','abstol',abstol,'alpha',alpha,
% 'n_sigma',n_sigma,fudge',fudge) |estimates the integral with integrand f
% to within an absolute error tolerance abstol with guaranteed uncertainty
% within alpha. All the field-value pairs are optional and can be supplied
% in different order.|
%
% Q = *cubMC_g*(f,interval,in_param) |estimates the integral with integrand f
% to within an absolute error tolerance in_param.abstol with guaranteed
% uncertainty within in_param.alpha. If a field is not specified, the
% default value is used.|
%
% *Input Arguments*
%
% * f --- |the integrand|
%
% * interval --- |the integration interval|
%
% * in_param.measure --- |the measure for generating the random variable,
%   the default is uniform.|
%
% * in_param.abstol --- |the absolute error tolerance, default value is
%   1e-2.|
%
% * in_param.alpha --- |the uncertainty, default value is 1%.|
%
% * in_param.n_sigma --- |initial sample size for estimating the sample
%   variance, the default value is 1e3.|
%
% * in_param.fudge --- |the standard deviation inflation factor, the
%                       default value is 1.1.|
%
% *Output Arguments*
%
% * Q --- |the estimated value of the the integration.|
%
% * out_param.time --- |the time eclipsed.|
%
% * out_param.exit --- |the state of program when exiting.|
%
%                          0  success
%
%                         10  interval does not contain numbers
%
%                         11  interval not 2 x d
%
%                         12  interval is only a point in one direction
%
%                         13  interval is infinite when measure is uniform
%
%                         14  interval is not doubly infinite when measure is normal
%
%
%
%
%
%% Examples
%   Example 1:
%
    f = @(x) x.^2;interval = [1;2]; Q = cubMC_g(f,interval,'abstol',1e-2)
%%
%   Example 2:
    f = @(x) x.^2;interval = [1;2]; Q = cubMC_g(f,interval)


%% See Also
%
% <html>  
% <a href="help_funappx_g.html">funappx_g</a>
% </html>
%
% <html>
% <a href="help_integral_g.html">integral_g</a>
% </html>
%
% <html>
% <a href="help_meanMC_g.html">meanMC_g</a>
% </html>
%
%% Reference
% [1]  N. Clancy, Y. Ding, C. Hamilton, F. J. Hickernell, and Y.~Zhang, 
%      The complexity of guaranteed automatic algorithms: Cones, not
%      balls, preprint, 2013, arXiv:1303.2412 [math.ST].