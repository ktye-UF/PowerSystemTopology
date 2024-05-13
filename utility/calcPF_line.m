function P = calcPF_line(mpc_this, line, V, A, B)
% [~, Ybus, ~, ~] = makeJac(mpc_this);
% B = full(imag(Ybus));
% % line 2-3
node1 = mpc_this.branch(line, 1);
node2 = mpc_this.branch(line, 2);
Vs = V(node1);
Vr = V(node2);
As = A(node1);
Ar = A(node2);
P = (Vs * Vr * B(node1,node2)) * sin(deg2rad(As) - deg2rad(Ar)) * mpc_this.baseMVA;
end