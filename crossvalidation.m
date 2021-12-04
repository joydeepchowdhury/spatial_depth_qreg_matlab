function optimum_h_or_neighborhood_size = crossvalidation(t_vector_X, X_static,...
    t_vector_Y, Y_static, method_for_h, type, Kernel)

% This function calculates the optimum bandwidth or neighborhood size using
% cross validation. If method_for_h = 1, it returns the optimum bandwidth,
% else it returns the optimum neighborhood size. X_static is the n-by-p
% data matrix of n covariate observations, Y_static is the n-by-q data
% matrix of the corresponding response observations. Three types of
% estimators can be used for cross validations, viz., the weighted
% pointwise mean, the weighted pointwise median and the weighted spatial
% median. To use the weighted pointwise mean, put type = 'pointwise_mean',
% to use the weighted pointwise median, put type = 'pointwise_median', and
% to use the weighted spatial median, put type = 'spatial_median'. The
% numbers n, p and q must be greater than 1. Kernel is a function handle
% representing the kernel function used to calculate the weights. 
% t_vector_X and t_vector_Y are the respective grids on which the 
% observations in X_static and Y_static are recorded.

sample_size = size(X_static,1);
p_covariate = 2;
p_response = 2;

X_distance = zeros(sample_size,sample_size);
for i=1:1:sample_size
    for j=1:1:sample_size
        if i < j
            t1 = X_static(i,:);
            t2 = X_static(j,:);
            if p_covariate < inf
                d = Lp_norm(t_vector_X, (t1 - t2), p_covariate);
            else
                d = max(abs(t1 - t2));
            end
            X_distance(i,j) = d;
        elseif i == j
            d = 0;
            X_distance(i,j) = d;
        else
            X_distance(i,j) = X_distance(j,i);
        end
    end
end

X_distance_upper_tringular_wo_diag = triu(X_distance,1);
X_distance_vector = reshape(X_distance_upper_tringular_wo_diag, 1, (sample_size * sample_size));
X_distance_vector_sorted = sort(X_distance_vector,'ascend');
X_distance_vector_sorted_positive = X_distance_vector_sorted(X_distance_vector_sorted > 0);
X_distance_limits = [X_distance_vector_sorted_positive(1), X_distance_vector_sorted_positive(end)];

if method_for_h == 1
    h_vector = linspace(X_distance_limits(1), X_distance_limits(2),100);
    h_vector_length = length(h_vector);
    
    h_vector_check = ones(size(h_vector));
    for i=1:h_vector_length
        h = h_vector(i);
        
        Indices_zero = (X_distance <= h);
        Indices_zero_row_sum = sum(Indices_zero,2);
        Indices_zero_row_sum_proper = Indices_zero_row_sum - 1;
        if sum(Indices_zero_row_sum_proper == 0) > 0
            h_vector_check(i) = 0;
        end
    end
    h_vector_proper = h_vector(h_vector_check > 0);
    h_vector_length_proper = length(h_vector_proper);
    
    Error_Type_temp_average = zeros(1,h_vector_length_proper);
    for i=1:h_vector_length_proper
        h = h_vector_proper(i);
        
        Type_temp = zeros(sample_size,size(Y_static,2));
        Error_Type_temp = zeros(1,sample_size);
        for j=1:1:sample_size
            Y = Y_static;
            X = X_static;
            distances = X_distance(j,:);
            distances(j) = [];
            target_Y = Y(j,:);
            target_X = X(j,:);
            Y(j,:) = [];
            X(j,:) = [];
            
            local_Y_values = Y((distances <= h),:);
            local_X_values = X((distances <= h),:);
            t_vector = 1:1:size(X,2);
            Weights = kernelweights(target_X, local_X_values, t_vector, h, Kernel);
            
            if strcmp(type, 'pointwise_median') == 1
                weighted_median = zeros(1,size(local_Y_values,2));
                for k=1:size(local_Y_values,2)
                    vector_concerned = local_Y_values(:,k);
                    [vector_concerned_sorted, vector_concerned_sorted_index] = ...
                        sort(vector_concerned,'ascend');
                    weights_by_sorted_index = Weights(vector_concerned_sorted_index);
                    cumulative_weights_by_sorted_index = cumsum(weights_by_sorted_index);
                    index_weighted_quantile = find(cumulative_weights_by_sorted_index >= 0.5, 1);
                    weighted_median(k) = vector_concerned_sorted(index_weighted_quantile);
                end
                Type_temp(j,:) = weighted_median;
            elseif strcmp(type, 'pointwise_mean') == 1
                local_Y_values_weighted = (Weights * ones(1,size(local_Y_values,2)))...
                    .* local_Y_values;
                Type_temp(j,:) = mean(local_Y_values_weighted,1);
            elseif strcmp(type, 'spatial_median') == 1
                Type_temp(j,:) = spatialquantile(local_Y_values, Weights, 0, 0, t_vector_Y);
            else
                error('error: Enter correct type.')
            end
            
            Error_Type_temp(j) = Lp_norm(t_vector_Y, (target_Y - Type_temp(j,:)), p_response);
        end
        Error_Type_temp_average(i) = mean(Error_Type_temp);
    end
    [~,optimum_h_index] = min(Error_Type_temp_average);
    optimum_h = h_vector_proper(optimum_h_index);
    
    optimum_h_or_neighborhood_size = optimum_h;
else
    Number_lower = 5;
    Portion_higher = 0.5;
    Number_higher = ceil(Portion_higher * sample_size);
    Neighborhood_size_vector = (Number_lower:1:Number_higher);
    
    Error_Type_temp_average = zeros(1,length(Neighborhood_size_vector));
    for i=1:length(Neighborhood_size_vector)
        neighborhood_size = Neighborhood_size_vector(i);
        
        Type_temp = zeros(sample_size,size(Y_static,2));
        Error_Type_temp = zeros(1,sample_size);
        for j=1:1:sample_size
            Y = Y_static;
            X = X_static;
            distances = X_distance(j,:);
            distances(j) = [];
            target_Y = Y(j,:);
            target_X = X(j,:);
            Y(j,:) = [];
            X(j,:) = [];
            
            sorted_Distance_X = sort(distances, 'ascend');
            h = sorted_Distance_X(neighborhood_size);
            
            local_Y_values = Y((distances <= h),:);
            local_X_values = X((distances <= h),:);
            t_vector = 1:1:size(X,2);
            Weights = kernelweights(target_X, local_X_values, t_vector, h, Kernel);
            
            if strcmp(type, 'pointwise_median') == 1
                weighted_median = zeros(1,size(local_Y_values,2));
                for k=1:size(local_Y_values,2)
                    vector_concerned = local_Y_values(:,k);
                    [vector_concerned_sorted, vector_concerned_sorted_index] = ...
                        sort(vector_concerned,'ascend');
                    weights_by_sorted_index = Weights(vector_concerned_sorted_index);
                    cumulative_weights_by_sorted_index = cumsum(weights_by_sorted_index);
                    index_weighted_quantile = find(cumulative_weights_by_sorted_index >= 0.5, 1);
                    weighted_median(k) = vector_concerned_sorted(index_weighted_quantile);
                end
                Type_temp(j,:) = weighted_median;
            elseif strcmp(type, 'pointwise_mean') == 1
                local_Y_values_weighted = (Weights * ones(1,size(local_Y_values,2)))...
                    .* local_Y_values;
                Type_temp(j,:) = mean(local_Y_values_weighted,1);
            elseif strcmp(type, 'spatial_median') == 1
                Type_temp(j,:) = spatialquantile(local_Y_values, Weights, 0, 0, t_vector_Y);
            else
                error('error: Enter correct type.')
            end
            
            Error_Type_temp(j) = Lp_norm(t_vector_Y, (target_Y - Type_temp(j,:)), p_response);
        end
        Error_Type_temp_average(i) = mean(Error_Type_temp);
    end
    [~,optimum_neighborhood_size_index] = min(Error_Type_temp_average);
    optimum_neighborhood_size = Neighborhood_size_vector(optimum_neighborhood_size_index);
    
    optimum_h_or_neighborhood_size = optimum_neighborhood_size;
end

end

function lp_norm = Lp_norm(t_vector_local, Data_for_norm, p_local)

lp_norm = (trapz(t_vector_local, abs(Data_for_norm).^p_local, 2)).^(1/p_local);

end