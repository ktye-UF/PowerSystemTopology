 % Copyright (C)  2019 Mengyang Gu
% 
% This file is a part of the RobustGaSP Package in Matlab
% 
% The R version of the RobustGaSP Package is available at CRAN
% 
% RobustGaSP is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 2 of the License, or
% (at your option) any later version.
% 
% Robus GaSP is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
%
%  Mengyang Gu, July 2019
%

function [f,g] = neg_lik_deriv_ppgasp(param,nugget, nugget_est, R0, X,zero_mean,output, kernel_type,alpha,CL,a,b)
% Calculate objective f

[f,g] = ppgasp_post(param,nugget, nugget_est, R0, X,zero_mean,output, kernel_type,alpha,CL,a,b);
f=-f;
if nargout > 1
g=-g.*exp(param);
end

%if nargout > 1 % gradient required
%    g = [-400*(x(2)-x(1)^2)*x(1)-2*(1-x(1));
%        200*(x(2)-x(1)^2)];
%end


