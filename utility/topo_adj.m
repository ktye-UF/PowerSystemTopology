function [X_adj, Y_adj] = topo_adj(dataset, opt)
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
    % % self-loop
    if opt.loop     % self-loop: diagonal = 1
        A  = A + eye(size(A,1));
    end
    % % embed method: none (direct matrix), eig, svd
    X_adj = [X_adj; topo_embed(x, A, opt)];
    Y_adj = [Y_adj; dataset.ctgc(i).Y];
end
