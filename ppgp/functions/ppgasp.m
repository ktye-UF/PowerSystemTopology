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
% RobustGaSP is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% construct a parallel partial Gaussian stochastic process emulator 
%
% Syntax:
%   model = ppgasp(design,response)
%   model = variogram(design,response,options)
%
% Description:
%   This function constructs the emulator of a vector of ouput variables or
%   a scalar output variable 
%
% Input:
%   design - a n x p matrix of input variables of the computer models. Each row
%            is vector of input with coordinates with p dimensions and assume one
%            has n runs of computer model outputs.
%   response - a n x k matrix of output variables of the computer models.
%            Each row is a scalar or vector output of the computer models. 
%
%   options is a structure variable that may contain the following fields:
%     trend - the mean/trend matrix of inputs. The default value is a  vector of ones. 
%     zero_mean - it has zero mean or not. The default value is false 
%                  meaning the mean is not zero. true value means the mean is zero.
%     nugget - numerical value of the nugget variance ratio. If nugget is 
%              equal to 0, it means there is either no nugget or the nugget is estimated. 
%              If the nugget is not equal to 0, it means a fixed nugget. 
%              The default value is 0. 
%     nugget_est - logical value. true value means nugget should be estimated
%                   and false value means nugget is fixed or not estimated. 
%                   The default value is false.
%     range_par - either empty or a vector. If it is empty, it means range 
%                  parameters are estimated; otherwise range parameters are
%                  given. The default value is empty.
%     prior_choice - the current version only supports the choice of
%                    'ref_approx', which uses the jointly robust prior. See
%                    the following literature.
%     a - prior parameter in the jointly robust prior. The default value is 0.2. 
%     b - prior parameter in the jointly robust prior. The default value is
%         {n^{-1/p}(a+p)} where n is the number of runs and p is the 
%         dimension of the input vector. 
%     kernel_type - a vector specifying the type of kernels of each 
%                   coordinate of the input. 'matern_3_2' and 'matern_5_2'
%                   are Matern correlation with roughness parameter 3/2 and
%                   5/2 respectively. 'pow_exp' is power exponential 
%                   correlation with roughness parameter alpha. If 'pow_exp'
%                   is to be used, one needs to specify its roughness parameter. 
%                   The default choice is 'matern_5_2'.
%     alpha - roughness parameters in the 'pow_exp' correlation functions. 
%             The default choice is a vector with each entry being 1.9.
%     lower_bound -   logical value. true value means the default lower 
%                     bounds of the inverse range parameters are used to 
%                     constrain the optimization and false value
%                     means the optimization is unconstrained.
%                     The default value is true and we also suggest to use
%                     false in various scenarios. 
%     max_eval - the maximum number of iterations to estimate the range and
%                nugget parameters.  
%     initial_values - a matrix of initial values of the kernel parameters
%                      to be optimized numerically, where each row of the 
%                      matrix contains a set of the log inverse range 
%                      parameters and the log nugget parameter. 
%     num_initial_values - the number of initial values of the kernel 
%                          parameters in optimization.
%
% Output:
%   model - structure variable with the following fields:
%     input - a n x p matrix of input variables of the computer models.
%     output - a  n x k matrix of output variables of the computer
%              models.
%     num_obs - the number of computer model runs (n). 
%     k - the number of coordinate of output vector in each computer model
%         run.
%     p - the number of coordinte of the input vector. 
%     X - a n x q matrix for the mean/trend function.
%     zero_mean - whether the emulator has zero-mean.
%     q - the dimension of the mean basis function. 
%     nugget - the nugget value.
%     nugget_est - logical value pf estimating the nugget or not.
%     range_par - the range parameter. If estimated, it should be the 
%                 inverse of the beta_hat parameter.
%     prior_choice - the current version support 'ref_approx' (the jointly
%                   robust prior)
%     a - prior parameter of the jointly robust prior.
%     b - prior parameter of the jointly robust prior. 
%     kernel_type - the type of kernels used for emulation. 
%     alpha - roughness parameters in the power exponential kernel.
%     lower_bound - logical value of having a lower bound for estimating
%                   the kernel parameters. 
%     max_eval - maximum evaluation when estimating the kernel parameters.
%     initial_values - a matrix where each row is a set of initial value in
%                      estimating the logarithm of the inverse range parameter.
%     num_initial_values - the number of initial values. 
%     kernel_type_num - the number of kernel (one for each coordinate in
%                       the input).
%     R0 - an array where the j-th matrix is
%          an absolute difference matrix of the j-th input vector
%     CL - a vector of the scale of the input used for the lower bound and the prior.
%     LB - a p x 1 vector of the lower bound for inverse range parameters beta.
%     log_post - the logorithm of the posterior (with regard to a
%                 normalizing constant).
%     beta_hat - estimated inverse range parameters.
%     L - a matrix with dimension n x n for the Cholesky decomposition of the
%         correlation matrix R, i.e. L* L'=R.
%     LX - a matrix with dimension q x q for the cholesky decomposition of the X'*inv(R)*X.
%     theta_hat - a matrix of estimated parameters in the mean functions. 
%     sigma2_hat - the estimated variance parameters. 
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

%ppgasp function
%I only implement the version with JR prior in this version

function model = ppgasp(design,response,varargin)



model=struct();

model.input=design;
model.output=response;
   
if(length(varargin)>1)
    error('please put the additional inputs into a struct variable');
end

%fprintf('length of varargin is %d \n',length(varargin));

if length(varargin)==1
   params=varargin{:};
else
    params=struct();
end


%fprintf('nugget_est is %s \n',params.nugget_est);

if(size(model.output,1)~= size(model.input,1))
    error('the number of rows in the input should be the same as the number of rows in the output');  
end

[model.num_obs,model.k]=size(model.output);
model.p=size(model.input,2);

if(~isfield(params,'trend'))
  model.X=ones(model.num_obs,1);
else
  model.X=params.trend;
end

if(~isfield(params,'zero_mean'))
  model.zero_mean=false;
else
  model.zero_mean=params.zero_mean;
end


if(model.zero_mean)
    model.q=0;
else
    model.q=size(model.X,2);
end

%check for nugget

if(~isfield(params,'nugget'))
    model.nugget=0;
else
    model.nugget=params.nugget;
end

if(~isfield(params,'nugget_est'))
    model.nugget_est=false;
else
    model.nugget_est=params.nugget_est;
end

%

if(~islogical(model.nugget_est))
    error('nugget.est should be boolean (either T or F)');
end

if(model.nugget~=0 & model.nugget_est==T)
    error('one cannot fix and estimate the nugget at the same time'); 
end
    
if(~isnumeric(model.nugget))
    error('nugget should be a numerical value');  
end


  
if(~isfield(params,'range_par'))
    model.range_par=[];
else
    model.range_par=params.range_par;
    if(model.nugget_est)
        error('The current version of the package does not support fixing range parameters while estimating the nugget');  
    end
    if((size(params.range_par,1)*size(params.range_par,2))~=size(model.input,2) )
      error('range parameters should have the same dimension of an input.');    
    end
end


%In this version, I only implement the approach using 
%%the jointly robust prior 

%if(isempty(param_all{8}))
    model.prior_choice="ref_approx";
%else
%    model.prior_choice=param_all{9};
%end

if(~isfield(params,'a'))
    model.a=0.2;
else
    model.a=params.a;
end

if(~isfield(params,'b'))
   model.b=1/model.num_obs^(1/model.p)*(model.a+model.p);
   %model.b=1/(model.num_obs*model.k)^(1/model.p)*(model.a+model.p);

else
   model.b=params.b;
end

if(~isfield(params,'kernel_type'))
   model.kernel_type="matern_5_2";
else
   model.kernel_type=params.kernel_type;
end

if(length(model.kernel_type)==1)
    model.kernel_type=repmat(model.kernel_type,model.p,1);
elseif(length(model.kernel_type)~=model.p)
    error('Please specify the correct number of kernels.');
end
         

if(~isfield(params,'alpha'))
   model.alpha=1.9*ones(model.p,1);
else
   model.alpha=params.alpha;
end

if(~isfield(params,'lower_bound'))
   model.lower_bound=true;
else
   model.lower_bound=params.lower_bound;
end

if(~isfield(params,'max_eval') )
   model.max_eval=max(30,20+5*model.p);
else
   model.max_eval=params.max_eval;
end

if(~isfield(params,'initial_values'))
   model.initial_values=[];
else
   model.initial_values=params.initial_values;
end

if(~isfield(params,'num_initial_values'))
   model.num_initial_values=2;
else
   model.num_initial_values=params.num_initial_values;
end

model.kernel_type_num=repmat(0,model.p,1);
for(i_p=1: model.p)
    if(model.kernel_type(i_p)=="matern_5_2")
       model.kernel_type_num(i_p)=3;
    elseif (model.kernel_type(i_p)=="matern_3_2")
      model.kernel_type_num(i_p)=2;
    elseif (model.kernel_type(i_p)=="pow_exp")
      model.kernel_type_num(i_p)=1;
    end
end

model.R0=zeros(model.num_obs,model.num_obs,model.p);

for(i_p=1: model.p)
  model.R0(:,:,i_p)=abs(model.input(:,i_p)-model.input(:,i_p)');
end

model.CL=zeros(model.p,1);
for(i_cl=1:model.p)
    model.CL(i_cl)=(max(model.input(:,i_cl))-min(model.input(:,i_cl)))/model.num_obs^(1/model.p);
end
%%finish checking in and verifying the input parameters

if(isempty(model.range_par))
    COND_NUM_UB=10^(16);
    
    
    %R0=model.R0;
    %p=model.p;
    %kernel_type_num=model.kernel_type_num
    %alpha=model.alpha;
    %nugget=model.nugget;
    
    fun = @(param) search_LB_prob(param, model.R0, COND_NUM_UB,model.p,...
    model.kernel_type_num,model.alpha,model.nugget);

    
    param_0=0;
    lb = -5;
    ub = 12;
    
    options = optimset('MaxFunEvals',100);

   [LB_min,fval] = fminbnd(fun,lb,ub,options);
   
   LB_prob=exp(LB_min)/(exp(LB_min)+1);
    
   LB = zeros(model.p,1);
    
    for( i_LB = 1:model.p)
        LB(i_LB)=log(-log(LB_prob)/(max(max(model.R0(:,:,i_LB)))));
    end

    if model.lower_bound
        if model.nugget_est
            model.LB=[LB;-Inf];
        else
            model.LB=LB;
        end
    else
        if model.nugget_est
            model.LB=repmat(-Inf,model.p+1,1);
        else
            model.LB=repmat(-Inf,model.p,1);
        end
    end
    
    fprintf('The upper bounds of the range parameters are ');
    fprintf([repmat('%f ', 1, size(model.LB', 2)) '\n'],1./exp(model.LB'));
    
    %%set up some initial starts
    if isempty(model.initial_values)
        beta_initial=zeros(model.num_initial_values,model.p);
        eta_initial=zeros(model.num_initial_values,1);
        %Nov 2021, thoughts: this may still be too small to start
        beta_initial(1,:)=50*exp(LB);
        eta_initial(1)=0.0001;
        if model.num_initial_values>1
            beta_initial(2,:)=(model.a+model.p)./(model.p.*model.CL.*model.b)/2;
            eta_initial(2)=0.0002;
        end
        
        if model.num_initial_values>2
            for(i_ini = 3:model.num_initial_values)
                rng(i_ini)
                beta_initial(i_ini,:)=10^3*rand(model.p,1)./model.CL;
                eta_initial(i_ini)=10^(-3)*rand(1);
            end
        end
        model.initial_values=[log(beta_initial) log(eta_initial)];
        model.log_post=-Inf;
        
        for(i_ini = 1:model.num_initial_values)
            %note here I make it a vector (p times 1 not 1 times p)
           if model.nugget_est
               ini_value=model.initial_values(i_ini,:)';
           else
               ini_value=model.initial_values(i_ini,1:model.p)';
           end
           
            fprintf('The initial values of range parameters are ');
            fprintf([repmat('%f ', 1,model.p) '\n'],1./exp(ini_value(1:model.p)'));
            fprintf('Start of the optimization %d: \n', i_ini);
            
            
            
            %%optimization
            options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
            options = optimoptions(options,'MaxIterations',model.max_eval);
            options = optimoptions(options,'StepTolerance',10^(-5));
            options = optimoptions(options,'Display','off');
            options = optimoptions(options,'HessianApproximation','lbfgs');
               
            
            fun = @(param) neg_lik_deriv_ppgasp(param,model.nugget, model.nugget_est,...
                model.R0, model.X,model.zero_mean,model.output,model.kernel_type_num,...
                model.alpha,model.CL,model.a,model.b);
            
            
            A = [];
            b = [];
            Aeq = [];
            beq = [];
            lb = model.LB;
            ub = [];
            nonlcon = [];
            
            %%add error and control, Nov 2021  

            ME=[];
            tic
            %%add, Nov 2021  
            try 
            [param_opt fval exitflag output_opt] = fmincon(fun,ini_value,A,b,Aeq,beq,lb,ub,nonlcon,options);
            catch ME
            end
            toc
            
            %%added Nov 2021
            if(isempty(ME)) %%no error
                if ~model.nugget_est
                    nugget_par=model.nugget;
                else
                    nugget_par=exp(param_opt(model.p+1));
                end

                fprintf('The number of interation is %d \n',output_opt.iterations);
                fprintf('The number of function evaluation is %d \n',output_opt.funcCount);
                fprintf('The log marginal posterior is %f \n',-fval);
                fprintf('The optimized range parameters are ');
                fprintf([repmat('%f ', 1,model.p) '\n'],1./exp(param_opt(1:model.p)'));
                fprintf('The optimized nugget parameter is %f \n',nugget_par);
                fprintf(output_opt.message);

                if (-fval)>model.log_post
                   model.log_post=-fval
                   if(model.nugget_est)
                      model.beta_hat = exp(param_opt(1:model.p));
                      model.nugget=nugget_par;
                   else
                      model.beta_hat = exp(param_opt(1:model.p));
                   end
                end
            end
        end
        % give the range parameter if estimated
        model.range_par=1./model.beta_hat;

    else       
        model.LB=rep(-Inf,model@p);
        model.beta_hat=1/model.range_par;
    end
    
   
    
     [model.L, model.LX, model.theta_hat, model.sigma2_hat]=construct_ppgasp(...
         model.beta_hat, model.nugget, model.R0, model.X, model.zero_mean,...
                              model.output,model.kernel_type_num,model.alpha); 
         



end
    
  


%function model = ppgasp(design,response,trend,zero_mean,nugget,
%                        nugget_est,range_par,a,b,kernel_type,alpha,
%                        lower_bound,max_eval,initial_values,num_initial_values)

%end


%varargin={[1 2], [],  b}