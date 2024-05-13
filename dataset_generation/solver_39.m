function varargout = solver_39(X, varargin)
% % % !TODO: varargin changable
% % power flow solver for miniwecc
% % Input: 
% X: n_sample * dim_input, [load, gen] ([loadp, loadq, genp])
% % Output:
% Y: n_sample * dim_output, [voltage magnitude, voltage angle]
% ctime_pf: CPU time
% is_converge: check if voltage between limits

%% determine basecase
current_dir = pwd;  % Get the current directory
[~, current_folder, ~] = fileparts(current_dir);  % Extract the folder name
%%%%%%%%%%%%%%%%%%%
% % version 1: basecase file
% casefile = 'case39.m';
% % % get initial mpc
% mpc = loadcase(casefile); 
% mpc_this = mpc;     % !TODO: not necessary?
% % version 2: load mpc
if strcmp(current_folder, 'dataset_generation')
        load('save/mpc_this')
    else
        load('dataset_generation/save/mpc_this')
    end
load('save/mpc_this')
%%%%%%%%%%%%%%%%%%
% % 
mpopt = mpoption('verbose', 0, 'out.all', 0);   % no print
% runpf('WECC240_HS_2018_Basecase_modified', mpopt)
% % determine load idx
data_load.idx = union(find(mpc_this.bus(:,3)), find(mpc_this.bus(:,4)));    % bus that has P or Q
data_load.ratio = mpc_this.bus(data_load.idx,4) ./ mpc_this.bus(data_load.idx,3);   % ratio = Q/P
% % determine gen idx
% % determine output
bus_output = 1:size(mpc_this.bus,1);
%% prepare input
% % get dim
dim_bus = size(mpc_this.bus,1);
[n, dim_input] = size(X);
dim_load = length(data_load.idx);
% [~, dim_gen] = size(data_gen.value);
% if dim_input ~= dim_load+dim_gen
%     error('Input dimension error: not equal to load + gen')
% end

%% run simulation one by one
V = NaN(n, size(mpc_this.bus,1));
A = NaN(n, size(mpc_this.bus,1));
branch_p = NaN(n, size(mpc_this.branch,1));
branch_q = NaN(n, size(mpc_this.branch,1));
is_converge = ones(n, 1);
% [X(i, 1:dim_load)', X(i, 1:dim_load)' .* data_load.load_ratio.value, X(i, 1:dim_load)' .* data_load.load_ratio.value./X(i, 1:dim_load)']
% [X(i, 1:dim_load)', X(i, dim_load+1:dim_load*2)', X(i, dim_load+1:dim_load*2)'./X(i, 1:dim_load)']
tic
for i = 1:n     % iter through sample
    % % version 1: Q as latent
%     % change load
%     mpc_this.bus(data_load.idx, 3) = X(i, 1:dim_load)'; % 3rd col: P
%     mpc_this.bus(data_load.idx, 4) = X(i, 1:dim_load)' .* data_load.load_ratio.value; % 4th col: Q = P * ratio
%     % change gen
%     mpc_this.gen(data_gen.idx, 2) = X(i, dim_load+1:dim_load+dim_gen)'; % 2nd col: P
    % % version 2: Q as input
    % change load
    mpc_this.bus(data_load.idx, 3) = X(i, 1:dim_load)'; % 3rd col: P
    mpc_this.bus(data_load.idx, 4) = X(i, dim_load+1:dim_load*2)'; % 4th col: Q as input
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % change load version 2: change all bus (zeros still zeros)
    mpc_this.bus(:, 3) = X(i, 1:dim_bus)'; % 3rd col: P
    mpc_this.bus(:, 4) = X(i, dim_bus+1:dim_bus*2)'; % 4th col: Q as input
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % change gen
%     mpc_this.gen(data_gen.idx, 2) = X(i, dim_load*2+1:dim_load*2+dim_gen)'; % 2nd col: P
    
    % get result
    result = runpf(mpc_this, mpopt);
    % extract voltage, power flow
    V(i,:) = result.bus(:,8)';     % 8th column: magnitude
    A(i,:) = result.bus(:,9)';     % 9th column: angle (degree)
    branch_p(i,:) = result.branch(:,14);   % 14th col: P
    branch_q(i,:) = result.branch(:,15);   % 14th col: Q
    % check convergence
    if any(V(i,:)<0.9) || any(V(i,:)>1.1)
        is_converge(i) = 0;
    end
end
ctime_pf = toc;
%% determine output
% % !TODO: not as expected
% % % volt magnitude
% for i=1:length(bus_output)
%     idx_volt(i,1) = find(mpc.bus(:,1)==bus_output(i));
% end
% % % volt angle
% idx_output = [idx_volt, idx_volt+dim_load];  % 
% Y_all = [V A];     % n_sample * dim_output
% Y = Y_all(:, idx_output);


% % 
if nargin>1
    output_type = varargin{1};
else
    output_type = 'voltage';
end
switch output_type
    case 'voltage'
        Y = V;
    case 'angle'
        Y = A;
    case 'pf'
%         Y = branch_p;
        Y = [V, A, branch_p];
    otherwise
        error('output type not found')
end
% Y = [V,A,branch_p,branch_q];
varargout = {Y, ctime_pf, is_converge};








