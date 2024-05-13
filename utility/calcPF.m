function P = calcPF(mpc_this, V, A)
% % get Ybus and
[~, Ybus, ~, ~] = makeJac(mpc_this);
B = full(imag(Ybus));
% % get branch (node1-node2)
branch = mpc_this.branch(:,[1,2]);
% % calculate line power flow
for i=1:size(V,1)
    for j=1:size(branch,1)
        P(i,j) = calcPF_line(mpc_this, j, V(i,:), A(i,:), B);
    end
end
        


















