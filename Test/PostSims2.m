
StepsNum = 1500;
for Si = 1:StepsNum
    ConfigName = frameName('config_####.mat',Si);
    ForcesName = frameName('forces_####.mat',Si);
    load(ConfigName)
    load(ForcesName)
    indexWall_1 = collision_indices(:,2) == -1;
    indexWall_2 = collision_indices(:,2) == -2;
    xyzPointsWall_1{Si} = x(collision_indices(indexWall_1,1),:);
    xyzPointsWall_2{Si} = x(collision_indices(indexWall_2,1),:);
    uvwWall_1{Si} = v(collision_indices(indexWall_1,1),:);
    uvwWall_2{Si} = v(collision_indices(indexWall_2,1),:);
end
for j = 1:StepsNum
    xzPointsWall_1{j} = [xyzPointsWall_1{j}(:,1) xyzPointsWall_1{j}(:,3)];
    xzPointsWall_2{j} = [xyzPointsWall_2{j}(:,1) xyzPointsWall_2{j}(:,3)];
    uwWall_1{j} = [uvwWall_1{j}(:,1) uvwWall_1{j}(:,3)];
    uwWall_2{j} = [uvwWall_2{j}(:,1) uvwWall_2{j}(:,3)];
end
    
[xMeshWall_1,zMeshWall_1,uMeshWall_1,wMeshWall_1,nMeshWall_1] = binMesh(xzPointsWall_1,uwWall_1,[-4 4 -0.5 7.5],1,80,0,0);
[xMeshWall_2,zMeshWall_2,uMeshWall_2,wMeshWall_2,nMeshWall_2] = binMesh(xzPointsWall_2,uwWall_2,[-4 4 -0.5 7.5],1,80,0,0);

figure
subplot(1,2,1)
plot(uMeshWall_1,zMeshWall_1,'bo-')
hold on
plot(uMeshWall_2,zMeshWall_2,'ro-')

subplot(1,2,2)
plot(wMeshWall_1,zMeshWall_1,'bo-')
hold on
plot(wMeshWall_2,zMeshWall_2,'ro-')

% Dimensionless
dzMeshWall = zMeshWall_1(end)-zMeshWall_1(end-1);
qWall_1 = sum(uMeshWall_1(~isnan(uMeshWall_1)))*dzMeshWall;
ubarWall_1 = sum(uMeshWall_1(~isnan(uMeshWall_1)))/sum(~isnan(uMeshWall_1));
hWall_1 = qWall_1/ubarWall_1;
% Theory
ieta = 0:0.05:1;
uTheory = 7/3-35/6*ieta.^(3/2)+ 7/2*ieta.^(5/2);

figure
plot(uMeshWall_1./ubarWall_1,(zMeshWall_1(end)-zMeshWall_1)./hWall_1,'bo-')
axis equal
axis ij
hold on
plot(uTheory,ieta,'k')

%% 
yMin = -2.2;
yMax = 2.2;t
yStep = 0.2;
TotalRange = yMin:yStep:yMax;
for Si = 1:999
    ConfigName = frameName('config_####.mat',Si);
    load(ConfigName)
    for j = 1:length(TotalRange)-1
        Range = [TotalRange(j) TotalRange(j+1)];
        indexK = x(:,2) > Range(1) & x(:,2) < Range(2);
        xzPoints{Si,j} = [x(indexK,1) x(indexK,3)];
        uw{Si,j} = [v(indexK,1) v(indexK,3)];
    end
end

for iLayer = 1:length(TotalRange)-1
    xzPoints_iLayer = xzPoints(:,iLayer)';
    uw_iLayer = uw(:,iLayer)';
    [xMeshiLayer,zMeshiLayer,uMeshiLayer,wMeshiLayer,nMeshiLayer] = binMesh(xzPoints_iLayer,uw_iLayer,[-4 4 -0.5 7.5],1,80,0,0);
    uMesh_iLayerArc{iLayer} = uMeshiLayer;
    xMesh_iLayerArc{iLayer} = xMeshiLayer;
    zMesh_iLayerArc{iLayer} = zMeshiLayer;
    wMesh_iLayerArc{iLayer} = wMeshiLayer;
    nMesh_iLayerArc{iLayer} = nMeshiLayer;
    yMesh_iLayerArc{iLayer} = 0.5*(TotalRange(iLayer)+TotalRange(iLayer+1))*ones(length(uMeshiLayer),1);
end
figure
plot(uMesh_iLayerArc{1},zMesh_iLayerArc{1})
hold on
for i = 2:length(uMesh_iLayerArc)
    plot(uMesh_iLayerArc{i},zMesh_iLayerArc{i})
end

uyzPoints = [];
for i = 1:length(TotalRange)-1
    uyzPoints = [uyzPoints; uMesh_iLayerArc{i} yMesh_iLayerArc{i} zMesh_iLayerArc{i}];
end
figure
plot3(uyzPoints(:,1),uyzPoints(:,2),uyzPoints(:,3),'b.')

[yMesh,zMesh] = meshgrid(-2.1:0.2:2.1,-0.5:0.2:7.5);
uMesh = griddata(uyzPoints(:,2),uyzPoints(:,3),uyzPoints(:,1),yMesh,zMesh);

figure
surf(uMesh,yMesh,zMesh)

uMesh_iLayerSum = zeros(length(uMeshiLayer),1);
for iLayer = 1:length(TotalRange)-1
    uMesh_iLayerSum = uMesh_iLayerSum + uMesh_iLayerArc{iLayer};
end
uMesh_iLayerAve= uMesh_iLayerSum./(length(TotalRange)-1);

figure
plot(uMeshWall_1,zMeshWall_1,'bo-')
hold on
plot(uMesh_iLayerAve,zMesh_iLayerArc{1},'ko-')