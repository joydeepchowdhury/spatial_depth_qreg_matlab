tic

load('Cigar')

state_code_vector = unique(Cigar(:,1));
counts_per_state = histc(Cigar(:,1),state_code_vector);

Number_states = length(state_code_vector);

Price_Cigar = zeros(Number_states,max(counts_per_state));
Population = zeros(Number_states,max(counts_per_state));
Population_above_16 = zeros(Number_states,max(counts_per_state));
CPI = zeros(Number_states,max(counts_per_state));
NDI = zeros(Number_states,max(counts_per_state));
Sales = zeros(Number_states,max(counts_per_state));
Min_Price = zeros(Number_states,max(counts_per_state));
for i=1:1:Number_states
    state_code = state_code_vector(i);
    state_block = Cigar((Cigar(:,1) == state_code),:);
    
    [~,indices] = sort(state_block(:,2),'ascend');
    state_block_sorted_by_year = state_block(indices,:);
    
    Price_Cigar(i,:) = state_block_sorted_by_year(:,3)';
    Population(i,:) = state_block_sorted_by_year(:,4)';
    Population_above_16(i,:) = state_block_sorted_by_year(:,5)';
    CPI(i,:) = state_block_sorted_by_year(:,6)';
    NDI(i,:) = state_block_sorted_by_year(:,7)';
    Sales(i,:) = state_block_sorted_by_year(:,8)';
    Min_Price(i,:) = state_block_sorted_by_year(:,9)';
end

toc

X_static = NDI;
Y_static = Sales;
coordinates_response = (1:1:size(Y_static,2)) + 62;
coordinates_covariate = (1:1:size(X_static,2)) + 62;
Response_name = 'Sales of packs of cigarette per capita';
Covariate_name = 'NDI per capita';
attached_string_response_plot = ' over years in the states';
attached_string_covariate_plot = ' over years in the states';

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
                d = (trapz(coordinates_covariate, (abs(t1 - t2)).^p_covariate, 2)).^(1/p_covariate);
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

toc

method_for_h = 1;
type = 'spatial_median'; % 'pointwise_mean','spatial_median','pointwise_median'
Kernel = @(z)ones(size(z));
manual = 0; % if manual = 0, then the bandwidth optimum_h is chosen automatically by cross-validation
if manual == 0
    optimum_h = crossvalidation(coordinates_covariate, X_static,...
        coordinates_response, Y_static, method_for_h, type, Kernel);
else
    optimum_h = 8000;
    optimum_neighborhood_size = 0.5 * sample_size;
end

toc

tic

L2_norm_X = (trapz(coordinates_covariate, (abs(X_static)).^2, 2)).^(1/2);
[~,indices_sorted_by_L2_norm_X] = sort(L2_norm_X,'ascend');

%% Computations for the slices and plots

numpics = 5;
Chosen_indices_prelim = ceil(linspace(3,sample_size-3, numpics));

Chosen_indices = indices_sorted_by_L2_norm_X(Chosen_indices_prelim);

Covariate_chosen = X_static(Chosen_indices,:);

Neighbourhood_size = zeros(numpics,1);
data_points_local = cell(numpics,1);
local_lower = zeros(numpics,size(Y_static,2));
local_upper = zeros(numpics,size(Y_static,2));

alpha = 0.05;
Spatial_Median = zeros(numpics,size(Y_static,2));
UpperBoundary = zeros(numpics,size(Y_static,2));
LowerBoundary = zeros(numpics,size(Y_static,2));
for i=1:1:numpics
    Y = Y_static;
    
    target = Chosen_indices(i);
    Distance_X = X_distance(target,:);
    target_Y = Y(target,:);
    
    if method_for_h == 1
        h = optimum_h;
    else
        sorted_Distance_X = sort(Distance_X, 'ascend');
        h = sorted_Distance_X(optimum_neighborhood_size);
    end
    
    W = (Distance_X <= h) / sum(Distance_X <= h);
    
    local_Y_values = Y(Distance_X <= h,:);
    data_points_local{i} = local_Y_values;
    
    Neighbourhood_size(i) = size(local_Y_values,1);
    
    if Neighbourhood_size(i) > 1
        local_lower(i,:) = min(local_Y_values);
        local_upper(i,:) = max(local_Y_values);
        
        local_Y_values_weights = ones(size(local_Y_values,1), 1);
        Spatial_Median(i,:) = spatialquantile(local_Y_values, local_Y_values_weights,...
            0, 0, coordinates_response);
        ConfidenceSet = spatialquantileconfidenceset(local_Y_values, local_Y_values_weights,...
            0, 0, coordinates_response, alpha);
        UpperBoundary(i,:) = ConfidenceSet(1,:) + Spatial_Median(i,:);
        LowerBoundary(i,:) = ConfidenceSet(2,:) + Spatial_Median(i,:);
    else
        local_lower(i,:) = local_Y_values;
        local_upper(i,:) = local_Y_values;
        
        Spatial_Median(i,:) = local_Y_values;
        UpperBoundary(i,:) = local_Y_values;
        LowerBoundary(i,:) = local_Y_values;
    end
end

toc

%% Figures for quantiles and depth regions

y1 = min([ min(min(Spatial_Median)), min(min(UpperBoundary)),...
    min(min(LowerBoundary)) ]);
y2 = max([ max(max(Spatial_Median)), max(max(UpperBoundary)),...
    max(max(LowerBoundary)) ]);
leeway = (y2 - y1) * 0.05;
y1 = y1 - leeway;
y2 = y2 + leeway;
y_limits = [y1, y2];
y_limits_spatial = y_limits;

y_limits_response = y_limits_spatial;

y1 = min(min(Covariate_chosen));
y2 = max(max(Covariate_chosen));
leeway = (y2 - y1) * 0.05;
y1 = y1 - leeway;
y2 = y2 + leeway;
y_limits = [y1, y2];
y_limits_covariates = y_limits;

figure
for i=1:1:2
    if i == 2
        for j=1:1:numpics
            subplot(2,numpics, ((i-1)*numpics + j))
            plot(coordinates_response,Spatial_Median(j,:),'LineWidth',2,'Color','k')
            hold all
            plot(coordinates_response,UpperBoundary(j,:),'--k')
            plot(coordinates_response,LowerBoundary(j,:),'--k')
            ylim(y_limits_response)
            xlim([min(coordinates_response), max(coordinates_response)])
            if j== 1
                ylabel({'Conditional Spatial Median','with Confidence Set'})
            end
            hold off
        end
    end
    if i == 1
        for j=1:1:numpics
            subplot(2,numpics, ((i-1)*numpics + j))
            plot(coordinates_covariate,Covariate_chosen(j,:),'k')
            ylim(y_limits_covariates)
            xlim([min(coordinates_covariate), max(coordinates_covariate)])
            if j== 1
                ylabel({'Covariate','Curves'})
            end
        end
    end
end
% figure
% for i=1:1:2
%     if i == 1
%         for j=1:1:numpics
%             subplot(2,numpics, ((i-1)*numpics + j))
%             plot(coordinates_covariate,Covariate_chosen(j,:),'k')
%             ylim(y_limits_covariates)
%             xlim([min(coordinates_covariate), max(coordinates_covariate)])
%             if j== 1
%                 ylabel({'Covariate','Curves'})
%             end
%         end
%     end
%     if i == 2
%         for j=1:1:numpics
%             subplot(2,numpics, ((i-1)*numpics + j))
%             plot(coordinates_response,Spatial_Median(j,:),'LineWidth',2,'Color','k')
%             hold all
%             for k=1:length(coordinates_response)
%                 step = 3 * (coordinates_response(end) - coordinates_response(1))...
%                     / length(coordinates_response);
%                 y = (LowerBoundary(j,k):1:UpperBoundary(j,k));
%                 x = coordinates_response(k) * ones(length(y), 1);
%                 plot(x, y, '.k')
%             end
%             ylim(y_limits_response)
%             xlim([min(coordinates_response), max(coordinates_response)])
%             if j== 1
%                 ylabel({'Conditional Spatial Median','with Confidence Set'})
%             end
%             hold off
%         end
%     end
% end