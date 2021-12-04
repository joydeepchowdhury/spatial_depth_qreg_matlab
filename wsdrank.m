function rankings = wsdrank(X_to_rank, X_data, X_data_weights, t_vector)

% This function calculates the ranks the observations in X_to_rank
% according to the decreasing value of the weighted spatial depth with
% respect to the data matrix X_data and corresponding weights in the column
% vector X_data_weights. The dimension of the data matrix X_data is n-by-p,
% where n is the sample size and p is the dimension of an observation
% vector. The dimension of the matrix X_to_rank is m-by-p, where m is the
% number of observations to rank. X_data_weights is a n-by-1 column vector,
% and t_vector is a row vector recording the grid on which the values of
% the functional observations underlying X_data and X_to_rank are recorded.
% The number p must be greater than 1.

number_of_points = size(X_to_rank,1);
wsd = zeros(number_of_points,1);
for i=1:number_of_points
    y = X_to_rank(i,:);
    
    difference = (ones(size(X_data,1),1) * y) - X_data;
    norm_difference = sqrt(trapz(t_vector, difference.^2, 2));
    check_nonzero_norm = (norm_difference ~= 0);
    
    if sum(check_nonzero_norm) == 0
        weighted_average = zeros(size(difference));
    else
        difference_proper = difference(check_nonzero_norm,:);
        norm_difference_proper = norm_difference(check_nonzero_norm,:);
        weights_proper = X_data_weights(check_nonzero_norm);
        scaled_difference_proper = difference_proper ./ ...
            ( norm_difference_proper * ones(1,size(difference_proper,2)) );
        scaled_difference_proper_weighted = ( weights_proper * ones(1,size(difference_proper,2)) )...
            .* scaled_difference_proper;
        weighted_average = sum(scaled_difference_proper_weighted,1) / sum(X_data_weights);
    end
    
    weighted_spatial_depth_y = 1 - sqrt(trapz(t_vector, weighted_average.^2));
    wsd(i) = weighted_spatial_depth_y;
end

[~,rankings] = sort(wsd, 'descend');
end