% tic
% 
% sample_size = 100;
% 
% coordinates_covariate = (0:0.01:1);
% resolution_covariate = length(coordinates_covariate);
% B = unifrnd(0,1, sample_size,1);
% Covariate = ( (B * ones(1,resolution_covariate))...
%     .* (ones(sample_size,1) * exp(coordinates_covariate)) );
% 
% p_covariate = 2;
% 
% if p_covariate < inf
%     Norm_covariate = (trapz(coordinates_covariate, (abs(Covariate)).^p_covariate, 2)).^(1/p_covariate);
% else
%     Norm_covariate = max(abs(Covariate),[], 2);
% end
% 
% coordinates_response = coordinates_covariate;
% resolution_respose = length(coordinates_response);
% stepsize_response = (1/(resolution_respose-1));
% Increments_response = normrnd(0, sqrt(stepsize_response), sample_size,(resolution_respose-1));
% Z = [zeros(sample_size,1), cumsum(Increments_response,2)];
% 
% hetero_or_skew = 1; % 1 for heteroscedasticity, 2 for skewness
% hetero_hom_type = 1; % type 1 for heteroscedastic, type 2 for homoscedastic
% if hetero_or_skew == 1
%     if hetero_hom_type == 1
%         Response = (Norm_covariate * ones(1,resolution_respose)) .* Z;
%     end
%     if hetero_hom_type == 2
%         c_for_Bt = (1/4);
%         Response = Covariate + (c_for_Bt * Z);
%     end
% end
% 
% toc

tic

hetero_hom_type = 1; % type 1 for heteroscedastic data, type 2 for homoscedastic data
if hetero_hom_type == 1
    load('qr1dataheter.mat')
end
if hetero_hom_type == 2
    load('qr1datahom.mat')
end

X_static = Covariate;
Y_static = Response;
Response_name = 'Response';
Covariate_name = 'Covariate';

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
    optimum_h = 0.5;
    optimum_neighborhood_size = 0.5 * sample_size;
end

toc

tic

L2_norm_X = (trapz(coordinates_covariate, (abs(X_static)).^2, 2)).^(1/2);
[~,indices_sorted_by_L2_norm_X] = sort(L2_norm_X,'ascend');

%% Computations for the spread measures and plots

numpics = sample_size;
Chosen_indices_prelim = ceil(linspace(1,sample_size, numpics));

Chosen_indices = indices_sorted_by_L2_norm_X(Chosen_indices_prelim);

boxplot_percentage = 50;
Boxplot_lower = zeros(numpics,size(Y_static,2));
Boxplot_upper = zeros(numpics,size(Y_static,2));
Neighbourhood_size = zeros(numpics,1);
number_data_points_boxplot = zeros(numpics,1);
diameter_spread = zeros(numpics,1);
diameter_spread_indices = zeros(numpics,2);
data_points_boxplot = cell(numpics,1);
data_points_local = cell(numpics,1);
Spatial_quantile_1_positive = zeros(numpics,size(Y_static,2));
Spatial_quantile_1_negative = zeros(numpics,size(Y_static,2));
Spatial_spread = zeros(numpics,1);
for i=1:1:numpics
    Y = Y_static;
    Y_indices = (1:1:sample_size);
    
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
    local_Y_indices = Y_indices(Distance_X <= h);
    data_points_local{i} = local_Y_values;
    
    Neighbourhood_size(i) = size(local_Y_values,1);
    
    if Neighbourhood_size(i) > 1
        number_data_points_boxplot(i) = ceil((Neighbourhood_size(i) * boxplot_percentage) / 100);
        local_Y_values_weights = ones(size(local_Y_values,1), 1);
        WSDrankings = wsdrank(local_Y_values, local_Y_values, local_Y_values_weights,...
            coordinates_response);
        local_Y_values_ranked = local_Y_values(WSDrankings,:);
        local_Y_values_ranked_indices = local_Y_indices(WSDrankings);
        local_Y_values_boxplot_data = local_Y_values_ranked(1:number_data_points_boxplot(i),:);
        local_Y_values_boxplot_indices = local_Y_values_ranked_indices(1:number_data_points_boxplot(i));
        Boxplot_lower(i,:) = min(local_Y_values_boxplot_data);
        Boxplot_upper(i,:) = max(local_Y_values_boxplot_data);
        data_points_boxplot{i} = local_Y_values_boxplot_data;
        
        pairwise_distance = zeros(number_data_points_boxplot(i),number_data_points_boxplot(i));
        pairwise_distance_indices = cell(number_data_points_boxplot(i),number_data_points_boxplot(i));
        for k=1:1:number_data_points_boxplot(i)
            for l=1:1:number_data_points_boxplot(i)
                if k < l
                    t1 = local_Y_values_boxplot_data(k,:);
                    t1_index = local_Y_values_boxplot_indices(k);
                    t2 = local_Y_values_boxplot_data(l,:);
                    t2_index = local_Y_values_boxplot_indices(l);
                    if p_response < inf
                        d = (trapz(coordinates_response, (abs(t1 - t2)).^p_response, 2)).^(1/p_response);
                    else
                        d = max(abs(t1 - t2));
                    end
                    pairwise_distance(k,l) = d;
                    pairwise_distance_indices{k,l} = [t1_index,t2_index];
                elseif k == l
                    d = 0;
                    pairwise_distance(k,l) = d;
                    t_index = local_Y_values_boxplot_indices(k);
                    pairwise_distance_indices{k,l} = [t_index,t_index];
                else
                    pairwise_distance(k,l) = pairwise_distance(l,k);
                    pairwise_distance_indices{k,l} = pairwise_distance_indices{l,k};
                end
            end
        end
        [temp1,indices_1] = max(pairwise_distance);
        [temp2,index_2] = max(temp1);
        location_pairwise_distance = [indices_1(index_2),index_2];
        diameter_spread(i) = max(max(pairwise_distance));
        temp3 = pairwise_distance_indices{indices_1(index_2),index_2};
        diameter_spread_indices(i,:) = temp3;
        
        c_quantile = 0.5;
        u_index = 1;
        Spatial_quantile_1_positive(i,:) = spatialquantile(local_Y_values, local_Y_values_weights,...
            u_index, c_quantile, coordinates_response);
        Spatial_quantile_1_negative(i,:) = spatialquantile(local_Y_values, local_Y_values_weights,...
            u_index, - c_quantile, coordinates_response);
        t1 = Spatial_quantile_1_positive(i,:);
        t2 = Spatial_quantile_1_negative(i,:);
        Spatial_spread(i) = (trapz(coordinates_response, (abs(t1 - t2)).^p_response, 2)).^(1/p_response);
    else
        Boxplot_lower(i,:) = local_Y_values;
        Boxplot_upper(i,:) = local_Y_values;
        data_points_boxplot{i} = local_Y_values;
        diameter_spread(i) = 0;
        
        Spatial_quantile_1_positive(i,:) = local_Y_values;
        Spatial_quantile_1_negative(i,:) = local_Y_values;
        Spatial_spread(i) = 0;
    end
end

Neighbourhood_size_all = Neighbourhood_size;

toc

sample_indices = (1:1:sample_size);

cutoff = 3;
sample_indices_adequate = sample_indices(Neighbourhood_size > cutoff);
diameter_spread_adequate = diameter_spread(Neighbourhood_size > cutoff);
Spatial_spread_adequate = Spatial_spread(Neighbourhood_size > cutoff);
figure
plot(sample_indices_adequate, diameter_spread_adequate, 'k',...
    sample_indices_adequate, Spatial_spread_adequate, 'k:')
l1 = legend('$$ \widehat{D}_1( p \,|\, \mathbf{x} ) $$', '$$ \widehat{D}_2( \tau \,|\, \mathbf{x} ) $$',...
    'Location','NorthWest');
set(l1, 'interpreter','latex')
ylim([0 2])
xlabel('Ranks of the covariate curves')
ylabel('Conditional Spread')

%% Computations for the slices and plots

numpics = 5;
Chosen_indices_prelim = ceil(linspace(1,sample_size, numpics));

Chosen_indices = indices_sorted_by_L2_norm_X(Chosen_indices_prelim);

Covariate_chosen = X_static(Chosen_indices,:);

boxplot_percentage = 50;
Boxplot_lower = zeros(numpics,size(Y_static,2));
Boxplot_upper = zeros(numpics,size(Y_static,2));
number_data_points_boxplot = zeros(numpics,1);
data_points_boxplot = cell(numpics,1);

Neighbourhood_size = zeros(numpics,1);
data_points_local = cell(numpics,1);
Pointwise_Mean_local = zeros(numpics,size(Y_static,2));
Pointwise_Median_local = zeros(numpics,size(Y_static,2));
local_lower = zeros(numpics,size(Y_static,2));
local_upper = zeros(numpics,size(Y_static,2));

Spatial_Median = zeros(numpics,size(Y_static,2));
Spatial_quantile_1_positive = zeros(numpics,size(Y_static,2));
Spatial_quantile_1_negative = zeros(numpics,size(Y_static,2));
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
        Pointwise_Mean_local(i,:) = mean(local_Y_values);
        Pointwise_Median_local(i,:) = median(local_Y_values);
        
        number_data_points_boxplot(i) = ceil((Neighbourhood_size(i) * boxplot_percentage) / 100);
        local_Y_values_weights = ones(size(local_Y_values,1), 1);
        WSDrankings = wsdrank(local_Y_values, local_Y_values, local_Y_values_weights,...
            coordinates_response);
        local_Y_values_ranked = local_Y_values(WSDrankings,:);
        local_Y_values_boxplot_data = local_Y_values_ranked(1:number_data_points_boxplot(i),:);
        Boxplot_lower(i,:) = min(local_Y_values_boxplot_data);
        Boxplot_upper(i,:) = max(local_Y_values_boxplot_data);
        data_points_boxplot{i} = local_Y_values_boxplot_data;
        
        Spatial_Median(i,:) = spatialquantile(local_Y_values, local_Y_values_weights, 0, 0,...
            coordinates_response);
        
        c_quantile = 0.5;
        u_index = 1;
        Spatial_quantile_1_positive(i,:) = spatialquantile(local_Y_values, local_Y_values_weights,...
            u_index, c_quantile, coordinates_response);
        Spatial_quantile_1_negative(i,:) = spatialquantile(local_Y_values, local_Y_values_weights,...
            u_index, - c_quantile, coordinates_response);
    else
        local_lower(i,:) = local_Y_values;
        local_upper(i,:) = local_Y_values;
        Pointwise_Mean_local(i,:) = local_Y_values;
        Pointwise_Median_local(i,:) = local_Y_values;
        
        Boxplot_lower(i,:) = local_Y_values;
        Boxplot_upper(i,:) = local_Y_values;
        data_points_boxplot{i} = local_Y_values;
        
        Spatial_Median(i,:) = local_Y_values;
        Spatial_quantile_1_positive(i,:) = local_Y_values;
        Spatial_quantile_1_negative(i,:) = local_Y_values;
    end
end

toc

%% Figures for local curves

y1 = min([ min(min(local_lower)), min(min(Pointwise_Median_local)),...
    min(min(Pointwise_Mean_local)), min(min(local_upper)) ]);
y2 = max([ max(max(local_lower)), max(max(Pointwise_Median_local)),...
    min(min(Pointwise_Mean_local)), max(max(local_upper)) ]);
leeway = (y2 - y1) * 0.05;
y1 = y1 - leeway;
y2 = y2 + leeway;
y_limits = [y1, y2];
y_limits_local = y_limits;

y1 = min(min(Covariate_chosen));
y2 = max(max(Covariate_chosen));
leeway = (y2 - y1) * 0.05;
y1 = y1 - leeway;
y2 = y2 + leeway;
y_limits = [y1, y2];
y_limits_covariates = y_limits;

figure
for i=1:1:4
    if i == 1
        for j=1:1:numpics
            subplot(4,numpics, ((i-1)*numpics + j))
            plot(coordinates_response,data_points_local{j}', 'Color',[0.5 0.5 0.5])
            ylim(y_limits_local)
            xlim([min(coordinates_response)-0.02, max(coordinates_response)])
        end
    end
    if i == 2
        for j=1:1:numpics
            subplot(4,numpics, ((i-1)*numpics + j))
            plot(coordinates_response,Pointwise_Mean_local(j,:),'Color','k')
            ylim(y_limits_local)
            xlim([min(coordinates_response)-0.02, max(coordinates_response)])
        end
    end
    if i == 3
        for j=1:1:numpics
            subplot(4,numpics, ((i-1)*numpics + j))
            plot(coordinates_response,Pointwise_Median_local(j,:),'Color','k')
            ylim(y_limits_local)
            xlim([min(coordinates_response)-0.02, max(coordinates_response)])
        end
    end
    if i == 4
        for j=1:1:numpics
            subplot(4,numpics, ((i-1)*numpics + j))
            plot(coordinates_covariate,Covariate_chosen(j,:),'k')
            ylim(y_limits_covariates)
            xlim([min(coordinates_covariate)-0.02, max(coordinates_covariate)])
        end
    end
end

%% Figures for quantiles and depth regions

y1 = min([ min(min(Boxplot_lower)), min(min(Boxplot_upper)) ]);
y2 = max([ max(max(Boxplot_lower)), max(max(Boxplot_upper)) ]);
leeway = (y2 - y1) * 0.05;
y1 = y1 - leeway;
y2 = y2 + leeway;
y_limits = [y1, y2];
y_limits_boxplot = y_limits;

y1 = min([ min(min(Spatial_Median)), min(min(Spatial_quantile_1_positive)),...
    min(min(Spatial_quantile_1_negative)) ]);
y2 = max([ max(max(Spatial_Median)), max(max(Spatial_quantile_1_positive)),...
    max(max(Spatial_quantile_1_negative)) ]);
leeway = (y2 - y1) * 0.05;
y1 = y1 - leeway;
y2 = y2 + leeway;
y_limits = [y1, y2];
y_limits_spatial = y_limits;

y1 = min(y_limits_boxplot(1), y_limits_spatial(1));
y2 = max(y_limits_boxplot(2), y_limits_spatial(2));
y_limits_response = [y1, y2];

y1 = min(min(Covariate_chosen));
y2 = max(max(Covariate_chosen));
leeway = (y2 - y1) * 0.05;
y1 = y1 - leeway;
y2 = y2 + leeway;
y_limits = [y1, y2];
y_limits_covariates = y_limits;

figure
for i=1:1:3
    if i == 2
        for j=1:1:numpics
            subplot(3,numpics, ((i-1)*numpics + j))
            plot(coordinates_response,Spatial_Median(j,:),'LineWidth',2,'Color','k')
            hold all
            plot(coordinates_response,Spatial_quantile_1_positive(j,:),'--k')
            plot(coordinates_response,Spatial_quantile_1_negative(j,:),':k')
            ylim(y_limits_response)
            xlim([min(coordinates_response)-0.02, max(coordinates_response)])
            if j== 1
                ylabel({'Regression','Quantiles'})
            end
            hold off
        end
    end
    if i == 3
        for j=1:1:numpics
            subplot(3,numpics, ((i-1)*numpics + j))
            plot(coordinates_response,data_points_boxplot{j}', 'Color',[0.5 0.5 0.5])
            ylim(y_limits_response)
            xlim([min(coordinates_response)-0.02, max(coordinates_response)])
            if j== 1
                ylabel({'Conditional','Max. Depth Sets'})
            end
        end
    end
    if i == 1
        for j=1:1:numpics
            subplot(3,numpics, ((i-1)*numpics + j))
            plot(coordinates_covariate,Covariate_chosen(j,:),'k')
            ylim(y_limits_covariates)
            xlim([min(coordinates_covariate)-0.02, max(coordinates_covariate)])
            if j== 1
                ylabel({'Covariate','Curves'})
            end
        end
    end
end