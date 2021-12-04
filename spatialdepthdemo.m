n = 20;
neighbourhood_size = 7;

normalMarkerSize = 40; % 36 default
specialMarkerSize = 40; % 70 taken

Data = normrnd(0,1, n,2);

SpatialDepth = zeros(n,1);
for i=1:n
    point = Data(i,:);
    Datatemp = Data;
    Datatemp(i,:) = [];
    
    t1 = ( ones(size(Datatemp,1), 1) * point ) - Datatemp;
    t2 = sqrt(sum(t1.^2, 2));
    SpatialDistr = mean( t1 ./ (t2 * ones(1,2)) );
    SpatialDepth(i) = 1 - sqrt(sum(SpatialDistr.^2));
end

[~,I] = sort(SpatialDepth, 'descend');
highestDepthPoint = Data(I(1),:);
highestDepthPointNeighbours = Data(I(2:(neighbourhood_size+1)),:);
lowestDepthPoint = Data(I(n),:);
lowestDepthPointNeighbours = Data(I(1:neighbourhood_size),:);

figure

subplot(1,2,1)

MarkerSize = normalMarkerSize * ones(n-1,1);
ScatterColor = ones(n-1,1) * [0 0 0];
% ScatterColor(I(1),:) = [0 0 0];
% MarkerSize(I(1),:) = specialMarkerSize;
scatter(Data(I(1),1),Data(I(1),2), specialMarkerSize, 'k', 'd', 'LineWidth',2)
hold on
Data_1 = Data;
Data_1(I(1),:) = [];
scatter(Data_1(:,1),Data_1(:,2), MarkerSize, ScatterColor, 'LineWidth',1)
hold on
for i=1:size(highestDepthPointNeighbours,1)
    head = highestDepthPointNeighbours(i,:);
    tail = highestDepthPoint;
    body = head - tail;
    quiver(tail(1),tail(2), body(1),body(2), 0, 'k', 'LineWidth',1)
end
text(Data(I(1),1)+0.07,Data(I(1),2)-0.12, 'A', 'HorizontalAlignment', 'center')
hold off

subplot(1,2,2)

MarkerSize = normalMarkerSize * ones(n-1,1);
ScatterColor = ones(n-1,1) * [0 0 0];
%ScatterColor(I(n),:) = [0 0 0];
%MarkerSize(I(n),:) = specialMarkerSize;
scatter(Data(I(n),1),Data(I(n),2), specialMarkerSize, 'k', 'd', 'LineWidth',2)
hold on
Data_2 = Data;
Data_2(I(n),:) = [];
scatter(Data_2(:,1),Data_2(:,2), MarkerSize, ScatterColor, 'LineWidth',1)
%hold on
for i=1:size(lowestDepthPointNeighbours,1)
    head = lowestDepthPointNeighbours(i,:);
    tail = lowestDepthPoint;
    body = head - tail;
    quiver(tail(1),tail(2), body(1),body(2), 0, 'k', 'LineWidth',1)
end
text(Data(I(n),1)-0.12,Data(I(n),2)+0.12, 'B', 'HorizontalAlignment', 'center')
hold off

% 
% p1 = [2 3];                         % First Point
% p2 = [9 8];                         % Second Point
% dp = p2-p1;                         % Difference
% 
% figure(1)
% quiver(p1(1),p1(2),dp(1),dp(2),0)
% grid
% axis([0  10    0  10])
% text(p1(1),p1(2), sprintf('(%.0f,%.0f)',p1))
% text(p2(1),p2(2), sprintf('(%.0f,%.0f)',p2))
