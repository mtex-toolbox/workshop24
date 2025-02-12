
% import the data
ebsd = EBSD.load("Y593_106A.ctf")

ebsd = ebsd.reduce

ebsd = ebsd.gridify;

%%

plot(ebsd('Aug'),ebsd('Aug').orientations,'figSize','huge')

%% Grain reconstruction using alpha shape

[grains, ebsd.grainId] = calcGrains(ebsd,'alpha',5)

grains = smooth(grains,5)

%% visualize the two main phases

figure(1)
plot(ebsd('Aug'),ebsd('Aug').orientations,'figSize','huge')
hold on
plot(grains('Aug').boundary,'lineWidth',2)
hold off

%%
figure(2)
plot(ebsd('Fo'),ebsd('Fo').orientations,'figSize','huge')
hold on
plot(grains('Fo').boundary,'lineWidth',2)
hold off

%% eliminate small grains

ebsd = ebsd(grains(grains.grainSize>5));

ebsd = ebsd.gridify


%% EBSD Denoising and Filling

F = halfQuadraticFilter
ebsdS = smooth(ebsd,F,'fill',grains)


%% visualize the two main phases

figure(3)
plot(ebsdS('Aug'),ebsdS('Aug').orientations,'figSize','huge')
hold on
plot(grains('Aug').boundary,'lineWidth',2)
hold off

%%
figure(4)
plot(ebsdS('Fo'),ebsdS('Fo').orientations,'figSize','huge')
hold on
plot(grains('Fo').boundary,'lineWidth',2)
hold off

%%

mori = grains.boundary('Aug','Aug').misorientation;

%%

clf
plotAngleDistribution(mori)

%%

mdf = calcDensity(mori,'halfwidth',5*degree)

%%

plot(mdf,'axisAngle')

%%

plotSection(mdf,'axisAngle',170*degree)

%%

annotate(tw)

%%

plotSection(mdf,'axisAngle')

%%


[value,tw] = max(mdf,'numLocal',2)

%%

axes = rotation.byAxisAngle(tw.CS.bAxis,linspace(0,180*degree,300)) * Miller(1,0,0,tw.CS);

moriLine = orientation.byAxisAngle(axes,180*degree)

plot(moriLine,'axisAngle')

%% line plot

clf
plot(mdf.eval(moriLine))

%%

angle(tw.axis('max'),round(tw.axis('max'))) / degree


%%

tw.axis

tw.angle ./ degree

%%

adf = calcAxisDistribution(mdf)

%%

plot(adf)

[value,modes] = max(adf,'numLocal',3)

annotate(modes)

%%

modes = orientation.byAxisAngle(modes,180*degree)

%%

plot(ebsdS('Aug'),ebsdS('Aug').orientations,'figSize','huge')
hold on
plot(grains('Aug').boundary,'lineWidth',2)
hold off


gB = grains.boundary('Aug','Aug');
isTwin1 = angle(gB.misorientation, tw(1)) < 6*degree;
isTwin2 = angle(gB.misorientation, tw(2)) < 6*degree;

hold on
plot(gB(isTwin1),'lineWidth',2,'linecolor','w')
plot(gB(isTwin2),'lineWidth',2,'linecolor','fuchsia')
hold off













