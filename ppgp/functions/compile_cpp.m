% Copyright (C)  2019 Mengyang Gu
% 
% This file is a part of the RobustGaSP Package available at Matlab
% 
% The R version of the RobustGaSP Package is available at CRAN:
% https://cran.r-project.org/web/packages/RobustGaSP/index.html
% 
% RobustGaSP is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 2 of the License, or
% (at your option) any later version.
% 
% RobustGaSP is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% compile the C++ files in the RobustGaSP matlab package.  
%
% Syntax:
%   compile_cpp();
%
%  Mengyang Gu, July 2019


function [] = compile_cpp()

path1 = ['-I' fullfile(pwd,'functions/Eigen')];
%path2 = ['-I' fullfile(pwd,'functions')];

%mex('-v',path1,path2,'ppgasp.cpp')  %%this  generate some verbal
%information
mex(path1,'functions/ppgasp_post.cpp');
mex(path1,'functions/separable_multi_kernel.cpp')
mex(path1,'functions/construct_ppgasp.cpp')
mex(path1,'functions/pred_ppgasp.cpp')

end