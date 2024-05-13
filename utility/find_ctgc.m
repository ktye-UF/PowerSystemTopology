function idx_ctgc = find_ctgc(mpc)
% % N-1 ctgc exclude those who will create island
% % [TODO]: include generation outage (currently manually operated)
n_branch = size(mpc.branch, 1);
idx_ctgc = [];
for i=1:n_branch
    % % outage
    mpc_this = mpc;
    mpc_this.branch(i,11) = 0;  % 0 means open
    % % get adjacency matrix
    edge = get_edge_from_mpc(mpc_this.branch);
    rows = edge(1,:)+1;     % edge_index follow python index, start from 0
    cols = edge(2,:)+1;
    vals = ones(1,size(edge,2));
    sparse_idx = {rows, cols, vals};
    A = full(sparse(sparse_idx{:}));
    % % check if there is more than 1 area
    % Convert the adjacency matrix to a graph object
    G = graph(A);
    % Compute the connected components of the graph object
    bins = conncomp(G);
    % Compute the number of areas
    num_areas = max(bins);
    % % Print the number of areas
    % disp(num_areas);
    if num_areas > 1
        idx_ctgc = [idx_ctgc, i];
    end
end












