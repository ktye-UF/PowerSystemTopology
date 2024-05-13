clearvars
addpath(genpath(pwd))
uqlab
current_dir = pwd;  % Get the current directory
[~, current_folder, ~] = fileparts(current_dir);  % Extract the folder name

%% parameters
% % parameter for input distribution
casefile = '39';
std_level = 0.15;
dist_type = 'Gaussian';
% % parameter for input data (X)
n_sample = 1000;
sample_method = 'lhs';
% % parameter for output
output_type = 'pf';
%%%%%%%%%%%%%%%%%%%%%%%
% % parameter for contingency
mpc = loadcase(['case', casefile]);
n_branch = size(mpc.branch, 1);
% % Note: some line outage will cause isolated node
idx_ctgc_iso = find_ctgc(mpc);
idx_ctgc = setdiff(1:n_branch, idx_ctgc_iso);
n_ctgc = length(idx_ctgc);


% % For 118
%% create input distribution
% % load info
% mpc = loadcase(['case', casefile]);
load_idx = union(find(mpc.bus(:,3)), find(mpc.bus(:,4)));    % bus that has P or Q
load_mean = mpc.bus(load_idx,3);   % P
load_ratio = mpc.bus(load_idx,4) ./ (mpc.bus(load_idx,3) + 1e-6); % ratio = Q/P
% % load parameter
dist_param.mean = load_mean;     % P
dist_param.std = std_level .* dist_param.mean;
dist_param.n_load = length(load_idx);

% %
myInput = initial_input_39(dist_type, dist_param);

%% generate input data for pf solver (Xpf)
% % get load value, P & Q, dim=21
X_loadp = uq_getSample(myInput, n_sample, sample_method);
X_loadq = X_loadp .* load_ratio';

% % get gen value
% % create X_bus: dim=39 (21 loads & the rest all 0)
X_busp = zeros(n_sample, size(mpc.bus,1));
X_busp(:, load_idx) = X_loadp;
X_busq = zeros(n_sample, size(mpc.bus,1));
X_busq(:, load_idx) = X_loadq;
X_bus = [X_busp, X_busq];

%% topology change: contingency (line outage)
% % Note: currently use same load value for all ctgc
% % modify btanch status
% mpc.branch(:,11) is branch status: 1 -> in-service; 0 --> out-of-service

dataset.idx_ctgc = idx_ctgc;
for i=1:length(idx_ctgc)
    % % change branch status
    mpc_this = mpc;
    mpc_this.branch(idx_ctgc(i),11) = 0;  % 0 means open
    if strcmp(current_folder, 'dataset_generation')
        save('save/mpc_this', 'mpc_this')
    else
        save('dataset_generation/save/mpc_this', 'mpc_this')
    end
    % % solve pf to get Y
    [Y, ctime_pf, is_converge] = solver_39(X_bus, output_type);
    % % store data
    dataset.ctgc(i).X = X_busp;
    dataset.ctgc(i).Y = Y;
    dataset.ctgc(i).idx_branch = idx_ctgc(i);
    dataset.ctgc(i).ctime_pf = ctime_pf;
    dataset.ctgc(i).is_converge = is_converge;
    % get & store edge data
    % Note: the edge_index is already in python index (start from 0)
    dataset.ctgc(i).edge_index = get_edge_from_mpc(mpc_this.branch);
end

%% save
error
if std_level == 0.1
    version = '_v1';
elseif std_level == 0.05
    version = '_v2';
elseif std_level == 0.15
    version = '_v4-2';
end
if strcmp(current_folder, 'dataset_generation')
    save(['save/data_ctgc_', casefile, version], 'dataset', 'myInput', 'mpc', 'dist_param', 'idx_ctgc')
else
    save(['dataset_generation/save/data_ctgc_', casefile, version], ...
        'dataset', 'myInput', 'mpc', 'dist_param', 'idx_ctgc')
end
% % v2: std=0.05; v2-2: Y=A
% % v3: std=0.02
% % v4: std=0.15, Y=[V,A,P]; v4-2: sampling = lhs
% % 

%%
var(Y)
histogram(Y(:,3))


