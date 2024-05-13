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
% fast prediction by a parallel partial Gaussian stochastic process emulator 
% for a calibration problem
%
% Syntax:
%   pred_model = predict_ppgasp(model,testing_input)
%   pred_model = predict_ppgasp(model,testing_input,options)
%
% Description:
%   This function makes prediction by the emulator of a vector of ouput variables or
%   a scalar output variable. It designs particularly for a calibration
%   when one needs to generate a large number runs of MCMC samples. A few
%   terms can be precomputed once rather than in every step of MCMC
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
% See also: predict_ppgasp
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


function pred_model = predict_ppgasp_matlab(model,testing_input,varargin)

% this term can be pre-computed if required
% For a calibration problem using MCMC for posterior samples, it only needs to be computed once

if model.zero_mean
    if ~isfield(model,'R_inv_output')
        opts1.LT = true;
        opts2.UT = true;
        model.R_inv_output=linsolve(model.L',linsolve(model.L,model.output,opts1),opts2);
    end
else
    if ~isfield(model,'R_inv_tilde_output') | ~isfield(model,'R_inv_X')   
        mean_est=model.X*model.theta_hat;
        tilde_output=model.output-mean_est;
        opts1.LT = true;
        opts2.UT = true;
        model.R_inv_tilde_output=linsolve(model.L',linsolve(model.L,tilde_output,opts1),opts2);
        model.R_inv_X=linsolve(model.L',linsolve(model.L,model.X,opts1),opts2);
    end
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


qt_025 = t_inv(0.025,(model.num_obs-model.q));
qt_975 = t_inv(0.975,(model.num_obs-model.q));

r0=zeros(num_testing,model.num_obs,model.p);

for(i_p=1:model.p)
  r0(:,:,i_p)=abs(testing_input(:,i_p)-model.input(:,i_p)');
end

r=separable_multi_kernel(r0,model.beta_hat,model.kernel_type_num,model.alpha);
     
     
if(model.zero_mean)
     pred_model.mean=r*model.R_inv_output;
     
       if(~mean_only)
            opts1.LT = true;
            opts2.UT = true;

           rt_R_inv= linsolve(model.L',linsolve(model.L,r',opts1),opts2)';
           c_star_star=zeros(num_testing,1);
           for( i_loc=1:num_testing)
              rtR_inv_r=rt_R_inv(i_loc,:)*r(i_loc,:)';
              c_star_star(i_loc)=1+model.nugget-rtR_inv_r;
           end

            pred_sigma_2_star=zeros(num_testing,model.k);


            for(loc_i=1:model.k)
              pred_sigma_2_star(:,loc_i)=  model.sigma2_hat(loc_i)*abs(c_star_star.array());
            end
            
            pred_model.lower95=pred_model.mean+sqrt(pred_sigma_2_star)*qt_025;
            pred_model.upper95=pred_model.mean+sqrt(pred_sigma_2_star)*qt_975;
            pred_model.sd=sqrt(pred_sigma_2_star*model.num_obs/(model.num_obs-2));
            %(pred_sigma_2_star*(num_obs)/(num_obs-2)).array().sqrt().matrix();
       end
else
       pred_model.mean=testing_trend*model.theta_hat+ r*model.R_inv_tilde_output;
       if(~mean_only)
           
          opts1.LT = true;
          opts2.UT = true;


          rt_R_inv= linsolve(model.L',linsolve(model.L,r',opts1),opts2)';
          c_star_star=zeros(num_testing,1);
          q=size(testing_trend,2); 
          
          for( i_loc=1:num_testing)
               
              X_testing_X_R_inv_r_i=testing_trend(i_loc,:)-r(i_loc,:)*model.R_inv_X;
              diff2=X_testing_X_R_inv_r_i*linsolve(model.LX',linsolve(model.LX,X_testing_X_R_inv_r_i',opts1),opts2)';
              
              rtR_inv_r=rt_R_inv(i_loc,:)*r(i_loc,:)';
              c_star_star(i_loc)=1+model.nugget-rtR_inv_r+diff2;
          end
          
          
          pred_sigma_2_star=zeros(num_testing,model.k);


          for(loc_i=1:model.k)
              pred_sigma_2_star(:,loc_i)=  model.sigma2_hat(loc_i)*abs(c_star_star);
          end
         
            pred_model.lower95=pred_model.mean+sqrt(pred_sigma_2_star)*qt_025;
            pred_model.upper95=pred_model.mean+sqrt(pred_sigma_2_star)*qt_975;
            pred_model.sd=sqrt(pred_sigma_2_star*model.num_obs/(model.num_obs-2));
       end
end

end
