addpath('functions');

% One only need to compile once. 
%compile_cpp();

% Example 1: emulating humanity data from the DIAMOND computer model where 
%            the output is the number of casualties over 5 days after catastrophes. 
addpath('data/humanity_data');

humanity_X = dlmread('humanity_X.txt',' ',1,0);
humanity_Y = dlmread('humanity_Y.txt',' ',1,0);

humanity_Xt = dlmread('humanity_Xt.txt',' ',1,0);
humanity_Yt = dlmread('humanity_Yt.txt',' ',1,0);

% parameters for the  ppgasp emulator. Here we assume an estimated nugget
% and intercept.
design=humanity_X;
response=humanity_Y;

options.nugget_est=true;
options.lower_bound=false;


model=ppgasp(design,response,options);



% regular prediction

pred_model=predict_ppgasp(model,humanity_Xt);

%%pre-computation for not zero-mean case
mean_est=model.X*model.theta_hat;

tilde_output=model.output-mean_est;

opts1.LT = true;

opts2.UT = true;

model.R_inv_tilde_output=linsolve(model.L',linsolve(model.L,tilde_output,opts1),opts2);

model.R_inv_X=linsolve(model.L',linsolve(model.L,model.X,opts1),opts2);

%%%new prediction
opt.mean_only=true;

pred_model2=predict_ppgasp_matlab(model,humanity_Xt,opt);

max(abs(pred_model2.mean-pred_model.mean))