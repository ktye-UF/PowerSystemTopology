function [X_adj, Y_adj] = topo_lap_eigmap(dataset, opt)
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
    %%%%%%%%%%%%%%%%%%
    A = sparse(sparse_idx{:});
    % Scale the adjacency matrix so the maximal distance is 1
    A = A.^2;
    max_distance = max(nonzeros(A));
    A = A./max_distance;
    % Build an undirected graph from the adjacency matrix
    G = graph(A);
    % Evaluate Gaussian kernel on nonzero elements of the adjacency matrix
    A = spfun(@(x) exp(-x / (2 * opt.std ^ 2)), A);
    % Construct the diagonal degree matrix
    D = diag(sum(A, 2));
    % Compute the (unnormalized) graph Laplacian
    L = D - A;
    L = full(L);

%     % % self-loop
%     if opt.loop     % self-loop: diagnal = 1
%         A  = A + eye(size(A,1));
%     end
%     % % get laplacian matrix
%     D = diag(sum(A));
%     L = D - A;
    % % normalize
    if opt.norm
        L = D^(-1/2) * L * D^(-1/2);
    end
    % % embed method: none (direct matrix), eig, svd
    X_adj = [X_adj; topo_embed(x, L, opt.embed)];
    Y_adj = [Y_adj; dataset.ctgc(i).Y];
end
