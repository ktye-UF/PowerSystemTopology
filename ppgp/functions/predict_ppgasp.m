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
% prediction by a parallel partial Gaussian stochastic process emulator 
% 
% Syntax:
%   pred_model = predict_ppgasp(model,testing_input)
%   pred_model = predict_ppgasp(model,testing_input,options)
%
% Description:
%   This function makes prediction by the emulator of a vector of ouput variables or
%   a scalar output variable 
%
% Input:
%   model - the structure variable produced by the ppgasp function.
%   testing_input - a matrix of input variables for predictions.
%
%   options is a structure variable that may contain the following fields:
%     testing_trend - the mean/trend function for testing. The number of 
%             mean basis should be the same in the constructed emulator.
%     mean_only - a logical value of producing the predictive mean only or
%                 with predictive credible interval.
%
% Output:
%   mean - the predictive mean.
%   lower95 - lower bound of the 95% posterior credible interval.
%   upper95 - upper bound of the 95% posterior credible interval.
%   sd - standard deviation of each testing_input.
%
%
% Example 1: emulating humanity data from the DIAMOND computer model where 
%            the output is the number of casualties over 5 days after catastrophes. 
%           
%     addpath('data/humanity_data');
%     addpath('functions');
% 
%     % One only need to compile once.
%     compile_cpp();
%     % Load the data. For each run of the computer model, the output has 5
%     %dimensions.
%
%     humanity_X = dlmread('humanity_X.txt',' ',1,0);
%     humanity_Y = dlmread('humanity_Y.txt',' ',1,0);
% 
%     humanity_Xt = dlmread('humanity_Xt.txt',' ',1,0);
%     humanity_Yt = dlmread('humanity_Yt.txt',' ',1,0);
% 
%     design=humanity_X;
%     response=humanity_Y;
%
%     % Here we assume an estimated nugget and intercept.
%
%     options.nugget_est=true;
%     options.lower_bound=false;
% 
%     model=ppgasp(design,response,options);
% 
%     % prediction
% 
%     pred_model=predict_ppgasp(model,humanity_Xt);
% 
% 
%     % one may trucate this example as they should be nonnegative
%     %pred_model.mean(find(pred_model.mean<0))=0;
% 
%     % predictive root of mean squared error
%     sqrt(mean((pred_model.mean-humanity_Yt).^2,'all'))
% 
%     % standard deviation
%     std(humanity_Yt,0,'all')
% 
%     % you may trucate predictive interval to be nonnegative for this example
%     % pred_model.lower95(find(pred_model.lower95<0))=0;
%     % pred_model.upper95(find(pred_model.upper95<0))=0;
% 
%     % number of held-out data covered in the nominal 95 precent predictive interval
%     size(find(humanity_Yt>=pred_model.lower95 ...
%         & humanity_Yt<=pred_model.upper95),1)/(size(humanity_Yt,1)*size(humanity_Yt,2))
%     % average length of confidence interval
%     mean(pred_model.upper95-pred_model.lower95,'all')
% 
%     % plot the hold out data and prediction with 95 precent predictive interval
% 
%     for i=1:5 
%         subplot_here=subplot(2,3,i);
% 
%         plot(1:120, humanity_Yt(:,i),'r.','DisplayName','test data')
%         ylabel('number of casualties')
%         xlabel('test number')
%         hold on 
%         legend(subplot_here)
% 
%         error_i=(pred_model.upper95(:,i)-pred_model.lower95(:,i))/2;
%         errorbar(1:120,pred_model.mean(:,i),error_i,'b.','DisplayName','prediction');
% 
%     end
%
% Example 2: a scalar output function with 2-dimensional input
% For function details and reference information, see:
% http://www.sfu.ca/~ssurjano/
%     nonpolynomial = @(x) ((30 + 5*x(:,1).*sin(5*x(:,1))).*(4 + exp(-5*x(:,2)))-100)/6;
%     n=50;
%     % We generate input by uniform deisng. A better design is the latin hypercube design. See lhsdesign().
%     x=[rand(n,1) rand(n,1)];
%     y=nonpolynomial(x);
% 
%     model=ppgasp(x,y);
% 
%     % test output
%     [x1_testing_mat x2_testing_mat]=meshgrid(0:0.02:1,0:0.02:1);
% 
%     num_testing=length(0:0.02:1)^2;
%     x_testing=zeros(num_testing,2);
%     x_testing(:,1)=reshape(x1_testing_mat,[num_testing,1]);
%     x_testing(:,2)=reshape(x2_testing_mat,[num_testing,1]);
% 
%     pred_model=predict_ppgasp(model,x_testing);
% 
%     y_testing=nonpolynomial(x_testing);
% 
%     % predictive root of mean squared error
%     sqrt(mean((pred_model.mean-y_testing).^2))
%     std(y_testing)
% 
%     size(find(y_testing>=pred_model.lower95 ...
%         & y_testing<=pred_model.upper95),1)/(size(y_testing,1)*size(y_testing,2))
%     % average length of confidence interval
%     mean(pred_model.upper95-pred_model.lower95,'all')
% 
%     % make a plot about test output and prediction 
%     y_testing_mat=reshape(y_testing, [sqrt(num_testing),sqrt(num_testing)]);
%     pred_model_mean_mat=reshape(pred_model.mean, [sqrt(num_testing),sqrt(num_testing)]);
% 
%     subplot(1,2,1);
% 
%     surf(x1_testing_mat,x2_testing_mat,y_testing_mat)
%     xlabel('input 1')
%     ylabel('input 2')
% 
%     view(2)
%     c1=colorbar;
%     ylabel(c1, 'test output')
%     caxis([0 10])
% 
% 
%     subplot(1,2,2);
% 
%     surf(x1_testing_mat,x2_testing_mat,pred_model_mean_mat)
%     xlabel('input 1')
%     ylabel('input 2')
%     view(2)
%     c2=colorbar;
%     ylabel(c2, 'prediction')
%     caxis([0 10])
%
%
% See also: ppgasp
%
% References:
%     M. Gu. and J.O. Berger (2016). Parallel partial Gaussian process emulation for computer models
%     with massive output. Annals of Applied Statistics, 10(3), 1317-1347.
%     M. Gu, X. Wang and J.O. Berger (2018), Robust Gaussian stochastic process emulation, Annals of
%     Statistics, 46(6A), 3038-3066.
%     M. Gu (2019), Jointly robust prior for Gaussian stochastic process in emulation, calibration and
%     variable selection, 14(3), Bayesian Analysis. 
%
% Mengyang Gu, July 2019

function pred_model = predict_ppgasp(model,testing_input,varargin)

if(length(varargin)>1)
    error('please put the additional inputs into a struct variable');
end

pred_model=struct();
if length(varargin)==1
   params=varargin{:};
else
    params=struct();
end

num_testing=size(testing_input,1);




if ~isfield(params,'testing_trend')
    testing_trend=ones(num_testing,1);
else
    testing_trend=params.testing_trend;
end

if ~isfield(params,'mean_only')
    mean_only=false;
else
    mean_only=params.mean_only;
end

        
if model.zero_mean
    pred_model.testing_trend=zeros(0,size(testing_input,1));
end

%Just in case you don't a stat tool box, I just a function online for
%percentile for t distribution
%qt_025 = tinv(0.025,(model.num_obs-model.q));
%qt_975 = tinv(0.975,(model.num_obs-model.q));

qt_025 = t_inv(0.025,(model.num_obs-model.q));
qt_975 = t_inv(0.975,(model.num_obs-model.q));

r0=zeros(num_testing,model.num_obs,model.p);

for(i_p=1:model.p)
  r0(:,:,i_p)=abs(testing_input(:,i_p)-model.input(:,i_p)');
end

  if(mean_only)
      
     pred_model.mean=pred_ppgasp(model.beta_hat,model.nugget,model.input,model.X,model.zero_mean,...
 model.output, testing_input,testing_trend,model.L,model.LX,model.theta_hat,...
                       model.sigma2_hat,qt_025,qt_975,r0,...
                           model.kernel_type_num,model.alpha,mean_only);
  else    
 [pred_model.mean, pred_model.lower95, pred_model.upper95, pred_model.sd]=...
     pred_ppgasp(model.beta_hat,model.nugget,model.input,model.X,model.zero_mean,...
 model.output, testing_input,testing_trend,model.L,model.LX,model.theta_hat,...
                       model.sigma2_hat,qt_025,qt_975,r0,...
                           model.kernel_type_num,model.alpha,mean_only);
  end
  
end

