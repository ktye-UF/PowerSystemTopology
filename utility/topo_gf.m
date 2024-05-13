function [X_adj, Y_adj] = topo_gf(dataset, opt)
n_ctgc = length(dataset.ctgc);
X_adj = [];
Y_adj = [];
for i=1:n_ctgc
    % % get x (original)
    x = dataset.ctgc(i).X;
    % % get adjacency matrix
    edge = dataset.ctgc(i).edge_index;
    rows = edge(1,:)+1;     % edge_index follow python index, start from 0
    cols = edge(2,:)+1;
    vals = ones(1,size(edge,2));
    sparse_idx = {rows, cols, vals};
    A = full(sparse(sparse_idx{:}));
    [U, V] = graph_factorization(A, opt);

    X_adj = [X_adj; x*U];
    Y_adj = [Y_adj; dataset.ctgc(i).Y];
end
