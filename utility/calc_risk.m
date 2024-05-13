% Assuming f and x are matrices, with each row corresponding to a set of data's KDE.
% Assuming the sizes of f and x are both m rows by n columns.
function probabilities = calc_risk(f, x, limit)
m = size(f, 1); % Number of sets of data
% limit = 1.5; % Define a threshold value
probabilities = zeros(m, 1); % Initialize a vector to store the probability for each set of data

% Loop through each set of data
for i = 1:m
    % Find the indices in the current row that exceed the limit
    over_limit_indices = x(i, :) > limit | x(i, :) < -limit;
    % Calculate the probability of exceeding the limit
    probabilities(i) = trapz(x(i, over_limit_indices), f(i, over_limit_indices));
end

% Display the results
% disp('Probabilities of exceeding the limit:');
% disp(probabilities);
