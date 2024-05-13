function expected_excesses = calc_excess(f, x, limit)
    m = size(f, 1); % Number of data sets
    expected_excesses = zeros(m, 1); % Initialize a vector to store the expected excess for each data set

    % Loop through each data set
    for i = 1:m
        % Find indices in the current row that exceed the limit
        over_limit_indices = x(i, :) > limit | x(i, :) < -limit;
        % Calculate the expected value for the part exceeding the limit
        if any(over_limit_indices)
            % Calculate the probability of exceeding the limit
            prob_excess = trapz(x(i, over_limit_indices), f(i, over_limit_indices));
            % Calculate the expected excess for the part exceeding the limit
            expected_excesses(i) = trapz(x(i, over_limit_indices), x(i, over_limit_indices) .* f(i, over_limit_indices)) / prob_excess;
        else
            expected_excesses(i) = 0;
        end
    end
    expected_excesses = abs(expected_excesses) - limit;
    % Display the results
    % disp('Expected excess for each data set:');
    % disp(expected_excesses);
end
