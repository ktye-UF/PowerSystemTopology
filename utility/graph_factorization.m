function [U, V] = graph_factorization(A, opt)
% 输入：
% A：邻接矩阵，大小为 n×n。
% k：嵌入空间的维度。
% lambda：正则化参数。
% max_iter：最大迭代次数。
% 输出：
% U：节点的向量表示，大小为 n×k。
% V：邻居节点的向量表示，大小为 n×k。

n = size(A, 1);
k = opt.k;
lambda = opt.lambda;
max_iter = opt.max_iter;

% 初始化节点的向量表示
U = rand(n, k);

% 对邻接矩阵进行归一化处理
D = diag(sum(A, 2));
L = D - A;
if opt.norm
    L = L ./ sqrt(n);
end
I = eye(n);

% 迭代更新节点向量表示
for iter = 1:max_iter
    V = (I + lambda * L) \ (lambda * A * U);
    U = (I + lambda * L) \ (lambda * A' * V);
end

end
