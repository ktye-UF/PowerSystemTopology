function edge = get_edge_from_mpc(Line)
Line(Line(:,11)==0,:) = [];     % remove out of service line
[m,n]=size(Line);
edge=[];
for i=1:m
    start_node=Line(i,1)-1;
    end_node=Line(i,2)-1;
    edge=[edge [start_node;end_node] [end_node; start_node]];
end