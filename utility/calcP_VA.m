function [P1, P2, mape, mape_br] = calcP_VA(X_p,Y_p,ctgc_this,model,mpc,idx_ctgc)
x_testing = X_p((ctgc_this-1)*1000+1:ctgc_this*1000,:);
% y_testing = Y_p((ctgc_this-1)*1000+1:ctgc_this*1000,idx_out);
y_testing_backup = Y_p((ctgc_this-1)*1000+1:ctgc_this*1000, :);

y_pred = predict_ppgasp_isotropic(model, x_testing);
ypred = y_pred.mean;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% together
%% from matpower
P0 = y_testing_backup(:,79:end);
%% param
% % change branch status
mpc_this = mpc;
mpc_this.branch(idx_ctgc(ctgc_this),11) = 0;  % 0 means open
%% test
V1 = y_testing_backup(:,1:39);
A1 = y_testing_backup(:,40:78);
P1 = calcPF(mpc_this, V1, A1);
%% predict
V = ypred(:,1:29);
V_gen = [1.0499,0.982,0.9841,0.9972,1.0123,1.0494,1.0636,1.0275,1.0265,1.03];
V_all = [V, repmat(V_gen, size(V,1), 1)];
A = ypred(:,30:58);
A_gen = y_testing_backup(:, 69:78);
A_all = [A, A_gen];
% % 
V2 = V_all;
A2 = A_all;
P2 = calcPF(mpc_this, V2, A2);

%% error
mape = abs((P1-P2)./(P1+1e-6)) * 100;
mape_br = mean(mape);
% mean(mean(mape))




















end