function Quantile = spatialquantile(Data_original, Weights, u_index, c, t_vector)

% This function computes the weighted spatial quantile from the data matrix
% X_data_original with corresponding weights in the variable Weights. The
% data matrix X_data_original is of dimension n-by-p and contains n 
% functional observations which are recorded on a grid of length p, 
% recorded in the row vector t_vector. The variable Weights is a column 
% vector of length n, whose entries are the weights of the corresponding
% observations in X_data_original. The variables u_index and c are a 
% positive integer and a real number in (-1, 1) respectively, which
% together determine the variable u for the spatial quantile computation.
% For example, if u_index = 1 and c = 0.5, then we compute the weighted
% spatial quantile corresponding to u equalling to 0.5 times the weighted 
% first principal component of the data in X_static_original.

%% Checking whether there is only one distinct observation in Data_original

z = Data_original(1,:);
Difference = ones(size(Data_original,1),1) * z - Data_original;
norm_Difference = sqrt(trapz(t_vector, Difference.^2, 2));
if sum(norm_Difference) == 0
    Quantile = Data_original;
    return
end

%% Dimension reduction for the input data

n = size(Data_original,1);
if sum(Weights) ~= 1
    Weights = Weights / sum(Weights);
end

t_1 = sqrt(n);
t_2 = 2 * n^(1/3);
t_3 = min(t_1,t_2);
d_n = floor(t_3);

Weighted_Mean = mean((Weights * ones(1,size(Data_original,2))) .* Data_original, 1);
Centred_Data = Data_original - ones(n,1) * Weighted_Mean;
Weighted_Cov_Matrix = Centred_Data' * diag(Weights) * Centred_Data;
Weighted_Cov_Matrix = (Weighted_Cov_Matrix + Weighted_Cov_Matrix') / 2;
[Eigenvectors, Eigenvalues] = eig(Weighted_Cov_Matrix);
vector_Eigenvalues = diag(Eigenvalues);
[~, index_Eigenvalues_sorted] = sort(vector_Eigenvalues,'descend');
Eigenvectors_sorted = Eigenvectors(:,index_Eigenvalues_sorted);
Coefficient_Matrix = Centred_Data * Eigenvectors_sorted;
Eigenvectors_sorted_truncated = Eigenvectors_sorted(:,1:d_n);
Coefficient_Matrix_truncated = Coefficient_Matrix(:,1:d_n);
Data_reduced = Coefficient_Matrix_truncated;

Data = Data_reduced;

if u_index == 0
    u = zeros(1,d_n);
elseif (u_index > 0) && (u_index <= d_n)
    temp = eye(d_n);
    u = c * temp(u_index,:);
else
    error('error: u_index must be smaller than d_n.')
end

%% Checking whether the weighted quantile is present in the data itself

Check = 0;
for i=1:n
    X = Data;
    x = X(i,:);
    U = X - ones(size(X,1),1) * x;
    weights_for_i = Weights;
    
    weighted_norms_U = weights_for_i .* sqrt(sum(U.^2,2));
    all_indices_for_i = 1:n;
    J_i = all_indices_for_i(weighted_norms_U == 0);
    J_i_complement = setdiff(all_indices_for_i, J_i);
    J_i = setdiff(J_i, i);
    
    U_new = U(J_i_complement,:);
    U_new = U_new ./ ( sqrt(sum(U_new.^2,2)) * ones(1,size(U_new,2)) );
    weights_proper = weights_for_i(J_i_complement);
    V = sum((weights_proper * ones(1,size(U_new,2))) .* U_new, 1) +...
        sum(weights_proper) * u;
    
    if sqrt(sum(V.^2)) <= (1 + sqrt(sum(u.^2))) * (sum(weights_for_i(J_i)))
        Quantile_coefficients = x;
        Check = 1;
        break
    end
end

%% Checking whether the data lie on a straight line, and computing the quantile in that case

x = Data(1,:);
for i=2:n
    y = Data(i,:);
    direction_vector = y - x;
    if sqrt(sum(direction_vector.^2)) > 0
        break
    end
end
Check_if_linear = zeros(1,n);
s = zeros(1,n);
for i=1:n
    z = Data(i,:);
    s_vector = (z - x) ./ direction_vector;
    if length(unique(s_vector)) == 1
        Check_if_linear(i) = 1;
        s(i) = unique(s_vector);
    end
end
if sum(Check_if_linear) == n
    projected_u = (u * direction_vector') / sqrt(sum(direction_vector.^2));
    alpha = (projected_u + 1) / 2;
    
    [s_sorted, s_sorted_index] = sort(s,'ascend');
    weights_sorted_index = Weights(s_sorted_index);
    cumulative_weights_sorted_index = cumsum(weights_sorted_index);
    index_weighted_quantile = find(cumulative_weights_sorted_index >= alpha, 1);
    s_weighted_quantile = s_sorted(index_weighted_quantile);
    
    Quantile_coefficients = x + (s_weighted_quantile * direction_vector);
    Check = 1;
end

%% Iteration procedure when the weighted quantile is not present in the data, or the data is not linear

if Check == 0
    X = Data;
    
    Q_1 = zeros(1,size(X,2));
    for i=1:length(u)
        vector_concerned = X(:,i);
        [vector_concerned_sorted, vector_concerned_sorted_index] = sort(vector_concerned,'ascend');
        weights_sorted_index = Weights(vector_concerned_sorted_index);
        cumulative_weights_sorted_index = cumsum(weights_sorted_index);
        index_weighted_quantile = find(cumulative_weights_sorted_index >= (u(i)+1)/2, 1);
        Q_1(i) = vector_concerned_sorted(index_weighted_quantile);
    end
    Q_best_till_now = Q_1;
    g_best_till_now = g_function_weighted(Data, Q_best_till_now, Weights, u);
    
    Phi = zeros(size(X,2), size(X,2));
    for i=1:n
        t1 = X(i,:) - Q_1;
        if sqrt(sum(t1.^2)) > 0
            Phi = Phi + Weights(i) * ...
                (( eye(size(X,2)) - ((t1' * t1) / sum(t1.^2)) ) / sqrt(sum(t1.^2)));
        end
    end
%     if cond(Phi) >= 10
%         warning('Bad initial value, final output may not be good estimate.')
%     end
    
    Threshold = 0.001;
    iteration_number = 1;
    maximum_iteration_number = 10000;
    while(iteration_number <= maximum_iteration_number)
        X = Data;
        Delta = zeros(1,size(X,2));
        Phi = zeros(size(X,2), size(X,2));
        for i=1:n
            t1 = X(i,:) - Q_1;
            if sqrt(sum(t1.^2)) > 0
                Delta = Delta + Weights(i) * (t1 / sqrt(sum(t1.^2)));
                Phi = Phi + Weights(i) * ...
                    (( eye(size(X,2)) - ((t1' * t1) / sum(t1.^2)) ) / sqrt(sum(t1.^2)));
            end
        end
        Delta = Delta + (sum(Weights) * u);
        Q_2 = Q_1 + (Phi \ (Delta'))';
        
        difference_relative = sqrt(sum((Q_2 - Q_1).^2)) /...
            ( max(sqrt(sum(Q_1.^2)), sqrt(sum(Q_2.^2))) );
        
        if difference_relative < Threshold
            Quantile_coefficients = Q_2;
            Check = 1;
            break
        else
            g_at_Q_1 = g_function_weighted(Data, Q_1, Weights, u);
            g_at_Q_2 = g_function_weighted(Data, Q_2, Weights, u);
            if g_at_Q_2 <= g_at_Q_1
                if g_at_Q_2 <= g_best_till_now
                    Q_best_till_now = Q_2;
                    g_best_till_now = g_function_weighted(Data, Q_best_till_now, Weights, u);
                end
            else
                Q_2 = (g_at_Q_2 * Q_1 + g_at_Q_1 * Q_2) / (g_at_Q_1 + g_at_Q_2);
                g_at_Q_2 = g_function_weighted(Data, Q_2, Weights, u);
                if g_at_Q_2 <= g_best_till_now
                    Q_best_till_now = Q_2;
                    g_best_till_now = g_function_weighted(Data, Q_best_till_now, Weights, u);
                end
            end
            
            iteration_counter = 1;
            while 1
                Phi_temp = zeros(size(X,2), size(X,2));
                for i=1:n
                    t2 = X(i,:) - Q_2;
                    if sqrt(sum(t2.^2)) > 0
                        Phi_temp = Phi_temp + Weights(i) * ...
                            (( eye(size(X,2)) - ((t2' * t2) / sum(t2.^2)) ) / sqrt(sum(t2.^2)));
                    end
                end
                
                if (cond(Phi_temp) <= 10) || (iteration_counter > 5)
                    break
                else
                    Q_2 = (Q_1 + Q_2) / 2;
                    iteration_counter = iteration_counter + 1;
                end
            end
            
            Q_1 = Q_2;
        end
        iteration_number = iteration_number + 1;
    end
end
if Check == 0
    Quantile_coefficients = Q_best_till_now;
end

%% Calculating the weighted quantile

Quantile = (Quantile_coefficients * Eigenvectors_sorted_truncated') + Weighted_Mean;

end

function g = g_function_weighted(X_local, Q_local, weights_local, u_local)

g = sum(weights_local .* ( sum(( (ones(size(X_local,1), 1) * Q_local) -...
    X_local ).^2, 2) ).^(1/2), 1) / sum(weights_local) - u_local * Q_local';

end