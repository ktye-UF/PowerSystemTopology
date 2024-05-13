function tsi = get_tsi(y)
% % y: N * T
for i=1:size(y,1)
%     theta = max(abs(y(i,:)-y(i,1)));
    theta = max(abs(y(i,:)));
    tsi(i) = 100 * (360-theta) / (360+theta);
end

