function [X_lap, Y_lap] = topo_lap(dataset, opt)
n_ctgc = length(dataset.ctgc);
X_lap = [];
Y_lap = [];
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
    if opt.loop     % self-loop: diagnal = 1
        A  = A + eye(size(A,1));
    end
    % % get laplacian matrix
    D = diag(sum(A));
    L = D - A;
    % % normalize
    if opt.norm
        L = D^(-1/2) * L * D^(-1/2);
    end
    % % embed method: none (direct matrix), eig, svd
    X_lap = [X_lap; topo_embed(x, L, opt)];
    Y_lap = [Y_lap; dataset.ctgc(i).Y];
end
