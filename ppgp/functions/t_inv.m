function y = t_inv(x,n);
% t_inv returns inverse cumulative function (quantile) of the student t distribution,
% which does not require the machine learning tool box in matlab
% 
% Syntax:
% 
% y = t_inv(p,v);
%
% Input:
%   p - the percentile.
%   v - the degrees of freedom.
%
% Output:
%   y - the quantile of the student t  distribution at p percentile with v
%       degrees of freedom.
%

% allocate output memory and check size of arguments
y = x+n-n;	
n = n+x-x;

y = (sign(x - 1/2).*sqrt(n./betaincinv(2*min(x, 1-x), n/2, 1/2) - n));
