% function [mean_exp, var_exp, std_exp, mae, rmae, mse] = get_exp_stat(Y)
% for i=1:length(Y)
%     mean_exp(i,:) = mean(Y{i});
%     var_exp(i,:) = var(Y{i});
%     std_exp(i,:) = std(Y{i});   
%     mae(i,:) = mean(abs(Y{1} - Y{i}));
%     rmae(i,:) = mean(abs((Y{1} - Y{i}) ./ Y{1}));
%     mse(i,:) = mean((Y{1} - Y{i}).^2);
% end
% end

function [error_mean_mae, error_var_mae] = get_exp_stat(Y)
for i=2:length(Y)
    mean_exp(i,:) = mean(Y{i});
    var_exp(i,:) = var(Y{i});
    std_exp(i,:) = std(Y{i});   
%     mae(i,:) = mean(abs(Y{1} - Y{i}));
%     rmae(i,:) = mean(abs((Y{1} - Y{i}) ./ Y{1}));
%     mse(i,:) = mean((Y{1} - Y{i}).^2);
    error_mean_mae(i-1,:) = abs(mean(Y{i}) - mean(Y{1}));
    error_var_mae(i-1,:) = abs(var(Y{i}) - var(Y{1}));
end
end