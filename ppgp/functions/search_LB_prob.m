% Copyright (C)  2019 Mengyang Gu
%
% This file is a part of the RobustGaSP Package in Matlab
% 
% The R version of the RobustGaSP Package is available at CRAN:
% https://cran.r-project.org/web/packages/RobustGaSP/index.html
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
% Search the default lower bound of range parameters. 
%
% Syntax:
%
% squared_diff=search_LB_prob(param, R0, COND_NUM_UB,p,kernel_type_num,alpha,nugget)
%
% Description: Function to find the values to construct the default lower bound of range parameters.
%
% Input: 
%  param: a vector of natural logarithm of inverse-range parameters and
%         natural logarithm of the nugget-variance ratio parameter.
%  R0: an array for an absolute difference matrix the input vector.
%  COND_NUM_UB: the maximum condition number of the correlation matrix.
%  p: the number of input parameters.
%  kernel_type: a vector specifying the type of kernels of each 
%               coordinate of the input. 'matern_3_2' and 'matern_5_2'
%               are Matern correlation with roughness parameter 3/2 and
%               5/2 respectively. 'pow_exp' is power exponential 
%               correlation with roughness parameter alpha. If 'pow_exp'
%               is to be used, one needs to specify its roughness parameter. 
%               The default choice is 'matern_5_2'.
%  alpha: roughness parameters in the power exponential kernel functions.
%  nugget: the nugget-variance ratio parameter if this parameter is fixed.
%
% Output:
%  squared_diff: a vector of values used in constructing the default lower 
%                bound of range parameters.
%
% Mengyang Gu, July 2019





function [squared_diff] = search_LB_prob(param, R0, COND_NUM_UB,p,kernel_type_num,alpha,nugget)
    num_obs=size(R0(:,:,1),1);
    propose_prob=exp(param)/(exp(param)+1);
    LB=zeros(p,1);
    for(i_LB=1:p)
        LB(i_LB)=log(-log(propose_prob)/(max(max(R0(:,:,i_LB)))));
    end
     R=separable_multi_kernel(R0,exp(LB),kernel_type_num,alpha);
     
     R=R+nugget*eye(num_obs);
     squared_diff=(cond(R)-COND_NUM_UB)^2;
end
