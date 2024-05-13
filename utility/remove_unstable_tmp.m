function [X, Y, idx] = remove_unstable_tmp(X, Y, type)
% % input Y: n * T
% % output y: same
x = X;
y = Y;
idx = [];
if size(y,2) >= 29
    num_out = 29;
else
    num_out = size(y,2);
end
switch type
    case 'angle'
        for i = 1:size(Y,1)
%             if any(Y(i,:)>1.2 | Y(i,:)<0.8 | any(isnan(Y(i,:))) )
            if any( any(isnan(Y(i,:))) )
                idx = [idx, i];
            end
        end
        x(idx,:) = [];
        y(idx,:) = [];
    case 'voltage'
        for i = 1:size(Y,1)
%             if any(Y(i,:)>1.2 | Y(i,:)<0.8 | any(isnan(Y(i,:))) )
            if any( Y(i,:)>1.1 | Y(i,:)<0.9 | any(isnan(Y(i,:))) )
                idx = [idx, i];
            end
        end
        x(idx,:) = [];
        y(idx,:) = [];
    case 'pf'
        for i = 1:size(Y,1)
%             if any(Y(i,:)>1.2 | Y(i,:)<0.8 | any(isnan(Y(i,:))) )
            if any( any(isnan(Y(i,:))) | Y(i,1:num_out)>1.2 | Y(i,1:num_out)<0.8 )
                idx = [idx, i];
            end
        end
        x(idx,:) = [];
        y(idx,:) = [];
    otherwise
        error('output type not found')
end
% % 
X = x;
Y = y;

