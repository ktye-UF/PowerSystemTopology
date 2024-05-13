clearvars
addpath(genpath(pwd))
% % load('dataset_ctgc')
%% parameter setting
output_type = 'pf';     % angle, voltage, power flow
% n_test = 10000;
n_batch = 20;       % number of sample for each ctgc/topology (50)      
% % 
switch output_type
    case 'angle'
        load('./dataset_generation/save/data_ctgc_39_v2-2')
        n_ctgc = length(dataset.ctgc);
        idx_gen = mpc.gen(:,1)';
        idx_out = 1:size(mpc.bus,1);     % has gen; [A]
    case 'voltage'
        load('./dataset_generation/save/data_ctgc_39_v2')
        n_ctgc = length(dataset.ctgc);
        idx_gen = mpc.gen(:,1)';
        idx_out = setdiff(1:size(mpc.bus,1), idx_gen);     % no gen; [V]
%         load('./dataset_generation/save/data_ctgc_118_v1')
    case 'pf'
        load('./dataset_generation/save/data_ctgc_39_v4-2')
        n_ctgc = length(dataset.ctgc);
        idx_gen = mpc.gen(:,1)';
        idx_out = [setdiff(1:size(mpc.bus,1), idx_gen), setdiff(1:39, idx_gen)+39];     % no gen; [V, A], no P
%         idx_out = 1:5;     % no gen; [V, A], no P
end

% idx_out = setdiff(1:size(mpc.bus,1), idx_gen);     % no gen

%% get X*A: A can be adjacency matrix, Laplacian matrix, eigenvector, etc.
matrix_type = 'lap';    % 'adj', 'lap': adjacency matrix, laplacian matrix
opt.embed = 'svd';      % 'none', 'eig', 'svd': none, eigenvector, svd
% % 
switch matrix_type
%%%%%%%%%%%%%%%%% scenario 1: adjacency matrix, eig or not
    case 'adj'
        opt.loop = true;        % self-loop -> little worse? not significant
%         opt.embed = 'none';     % 'none', 'eig', 'svd' 
        % get new input X_p
        [X_p, Y_p] = topo_adj(dataset, opt);

%%%%%%%%%%%%%%%%%%% scenario 2: laplacian matrix
    case 'lap'
        opt.loop = false;       % self-loop
        opt.embed = 'svd';    
        opt.norm = false;       % normalization
        opt.trunc = false;       % truncation (svd)
        [X_p, Y_p] = topo_lap(dataset, opt);

%%%%%%%%%%%%%%% scenario 3: laplacian eigenmap
% opt.loop = false;     
% opt.embed = 'svd';  
% opt.norm = false; 
% opt.std = 0.2;
% [X_p, Y_p] = topo_lap_eigmap(dataset, opt);

%%%%%%%%%%%%%%%%%%%% scenario 4: graph factorization
% % % not right gf
% opt.norm = true; 
% opt.k = 30;
% opt.lambda = 0.1;
% opt.max_iter = 1;
% [X_p, Y_p] = topo_gf(dataset, opt);

    otherwise 
        error_model('matrix type not found')
end

%% prepare input & output, traing & test
sample_method = 3;
switch sample_method 
% % version 1: sample n_batch for each topology  
    case 1      
        x = [];
        y = [];
        idx_tr = [];
%         for i=1:30    % transfer?
        for i=1:n_ctgc
            idx_batch = datasample((i-1)*1000+1:i*1000, n_batch, 'Replace', false);
            x = [x; X_p(idx_batch, :)];
            y = [y; Y_p(idx_batch, idx_out)];
            idx_tr = [idx_tr, idx_batch];
        end
        % % remove nan
        [x, y, idx_unstable] = remove_unstable_tmp(x, y, output_type);
        % % test samples
        idx_notr = setdiff(1:size(X_p,1), idx_tr);
        idx_test = datasample(idx_notr, n_test);
        x_testing = X_p(idx_test,:);
        y_testing = Y_p(idx_test,idx_out);
        % % remove nan
        [x_testing, y_testing, idx_unstable] = remove_unstable_tmp(x_testing, y_testing, output_type);
% % version 2: random
    case 2
        % % train
        [x, idx_tr] = datasample(X_p, n_batch*n_ctgc, 'Replace', false);
        y = Y_p(idx_tr, idx_out);
        y_backup = Y_p(idx_tr, :);  % used to get corresponding power flow, including 'idx_unstable_tr'
        [x, y, idx_unstable_tr] = remove_unstable_tmp(x, y, output_type);
        idx_notr = setdiff(1:size(X_p,1), idx_tr);
        idx_test = datasample(idx_notr, n_test);
        % % test
        x_testing = X_p(idx_test, :);
        y_testing = Y_p(idx_test, idx_out);
        y_testing_backup = Y_p(idx_test, :);    % used to get corresponding power flow
        [x_testing, y_testing, idx_unstable_te] = remove_unstable_tmp(x_testing, y_testing, output_type);
% % version 3: sample n_batch for each topology (same for test data)
    case 3
        x = [];
        y = [];
        % % train
        for i=1:n_ctgc
            idx_batch = (i-1)*1000+1:i*1000;
            x = [x; X_p(idx_batch(1:n_batch), :)];
            y = [y; Y_p(idx_batch(1:n_batch), idx_out)];
        end
        [x, y, idx_unstable_tr] = remove_unstable_tmp(x, y, output_type);
        % % test
        for i=1:n_ctgc
            idx_batch = (i-1)*1000+1:i*1000;
            x_testing = [x; X_p(idx_batch(201:1000), :)];
            y_testing = [y; Y_p(idx_batch(201:1000), idx_out)];
        end
        [x_testing, y_testing, idx_unstable_te] = remove_unstable_tmp(x_testing, y_testing, output_type);
    otherwise
        error_model('sample method not found')
end


%% PPGP construction
% % train
options.nugget_est = true;
options.lower_bound = false;
% options.isotropic = true;

tic
model = ppgasp(x,y,options);
% model = ppgasp_isotropic(x,y,options);
ctime_train = toc

% % predict
tic
y_pred = predict_ppgasp(model, x_testing);
% y_pred = predict_ppgasp_isotropic(model, x_testing);
ypred = y_pred.mean;
ctime_test = toc

% % error
mape = abs((ypred - y_testing) ./ y_testing);
error_model = abs(ypred - y_testing);

% mean(mean(mape))
% mean(error_model)
mean(mean(error_model))

% sqrt(sum((ypred - y_testing).^2))/size(y_testing,1)
error('-------save------')
%% save
% save(['save/ppgp_topo_v4-3_', matrix_type, opt.embed, '_nb', num2str(n_batch)])

%%%%%%%%%%%%%%%% End of Code %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Analysis %%%%%%%%%%%%%%%%%%%%%%%%
% % check function for usage
% calc_risk
% calc_excess
%%%%%%%%%%%%%%%% Plot figres just for illustration %%%%%%%%%%
%% plot all ctgc (selected branch output)
idx_plot = 1:10;
close all
for i=idx_plot
    figure; hold on;
    [f(1,:), xi(1,:)] = ksdensity(y_testing(:,i));
    plot(xi(1,:), f(1,:));
    [f(2,:), xi(2,:)] = ksdensity(ypred(:,i));
    plot(xi(2,:), f(2,:));
    legend('MC', 'PPGP')
    xlabel('Voltage magnitude/angle (pu/deg)'); ylabel('Probability density');
    title(['PDF of voltage magnitude/angle at bus '])
end

%% plot each bus different topology (all data)
ctgc_all = [1,2,3,4];
idx_plot = 1;
close all
for i=1:length(ctgc_all)
    x_testing = X_p((ctgc_all(i)-1)*1000+1:ctgc_all(i)*1000,:);
    y_testing = Y_p((ctgc_all(i)-1)*1000+1:ctgc_all(i)*1000,idx_out);
    [x_testing, y_testing, idx_unstable] = remove_unstable_tmp(x_testing, y_testing, output_type);
    y_pred = predict_ppgasp(model, x_testing);
    ypred = y_pred.mean;
    error_model = abs(ypred - y_testing);
    mean_error(i,:) = mean(error_model);
    figure; hold on;
    [f(1,:), xi(1,:)] = ksdensity(y_testing(:,idx_plot));
    plot(xi(1,:), f(1,:));
    [f(2,:), xi(2,:)] = ksdensity(ypred(:,idx_plot));
    plot(xi(2,:), f(2,:));
    legend('MC', 'PPGP')
    xlabel('Voltage magnitude/angle (pu/deg)'); ylabel('Probability density');
    title(['PDF of voltage magnitude/angle at bus ', num2str(idx_plot), ' with ctgc = ', num2str(ctgc_all(i))]);
end






















