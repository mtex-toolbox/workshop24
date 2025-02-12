%%                        MTEX Workshop 2024
%
%%   Exercise 10 - Orientation Relationships and parent grain reconstruction


%% Import some meteorite data

plotx2east

% import the ebsd data
mtexdata emsland

% extract crystal symmetries
cs_bcc = ebsd('Fe').CS;
cs_aus = ebsd('Aus').CS;

% recover grains
ebsd = ebsd('indexed');

[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',5*degree);
ebsd(grains(grains.grainSize<=2)) = [];
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',5*degree);

grains = smooth(grains,4);

%% visualize the data


plot(ebsd('Fe'),ebsd('Fe').orientations)
hold on
plot(grains.boundary,'lineWidth',2,'lineColor','gray')
plot(grains('Aus'),'FaceColor','blue','edgeColor','b','lineWidth',1,'DisplayName','Austenite')
hold off


%% plot the austenite orientations in axis angle


%% determine the mean austenite orientation

parentOri = mean(ebsd('Aus').orientations)

%% plot Ferrite orientations in pole figures

childOri = grains('Fe').meanOrientation;

h_bcc = Miller({1,0,0},{1,1,0},{1,1,1},cs_bcc);
h_fcc = Miller({1,0,0},{1,1,0},{1,1,1},cs_aus);

plotPDF(childOri,h_bcc,'MarkerSize',5,'MarkerFaceAlpha',0.05,'MarkerEdgeAlpha',0.1,'points',500);

nextAxis(1)
hold on
plot(parentOri * h_fcc(1).symmetrise ,'MarkerFaceColor','r')
xlabel('\((100)\)','Color','red','Interpreter','latex')

nextAxis(2)
plot(parentOri * h_fcc(3).symmetrise ,'MarkerFaceColor','r')
xlabel('\((111)\)','Color','red','Interpreter','latex')

nextAxis(3)
plot(parentOri * h_fcc(2).symmetrise ,'MarkerFaceColor','r')
xlabel('\((110)\)','Color','red','Interpreter','latex')
hold off

%%

% Kurdjumov Sachs



p2c = orientation.map(Miller(1,1,1,cs_aus),Miller(0,1,1,cs_bcc),...
      Miller(-1,0,1,cs_aus),Miller(-1,-1,1,cs_bcc))

childOri2 = parentOri * inv(p2c)

plotPDF(childOri2,h_bcc,'MarkerFaceColor','none','MarkerEdgeColor','k','linewidth',2)


%%

KS = orientation.KurdjumovSachs(cs_aus,cs_bcc)

%%

p2c = orientation.map(Miller(1,1,1,cs_aus),Miller(0,1,1,cs_bcc),...
      Miller(-1,0,1,cs_aus),Miller(-1,-1,1,cs_bcc))

plotPDF(variants(p2c,parentOri),'add2all','MarkerFaceColor','none','MarkerEdgeColor','k','linewidth',2)

%% How good are we with the Kurdjumov Sachs OR?

% Each parent-to-child misorientation can be calculated by
mori = inv(childOri) * parentOri;



%% Estimate the parent to child orientation relation ship

% The mean of all measured parent-to-child misorientations
p2cMean = mean(mori,'robust')

plotPDF(childOri,h_bcc,'MarkerSize',5,'MarkerFaceAlpha',0.05,'MarkerEdgeAlpha',0.1,'points',500);
hold on
plotPDF(variants(p2cMean,parentOri),'add2all','MarkerFaceColor','none','MarkerEdgeColor','k','linewidth',2)
hold off

% mean angular deviation in degree
mean(angle(mori, p2cMean)) ./ degree

%% Estimate only from child orientations

% extract all child to child misorientations
grainPairs = neighbors(grains('Fe'));
ori = grains(grainPairs).meanOrientation;

p2cIter = calcParent2Child(ori,KS)

%%


plot(childOri,'axisAngle')


%% Classification of the variants


[variantId, packetId] = calcVariantId(parentOri,childOri,p2cIter,'morito');

% colorize the orientations according to the variantID
color = ind2color(packetId);

plotPDF(childOri,color,h_bcc,'MarkerSize',5);


%%

plot(grains('Fe'),color)


%% Parent grain reconstruction
%% ---------------------------------------------------------------------


% load the data
mtexdata martensite
plotx2east

% grain reconstruction
[grains,ebsd.grainId] = calcGrains(ebsd('indexed'), 'angle', 3*degree);

% remove small grains
ebsd(grains(grains.grainSize < 3)) = [];

% reidentify grains with small grains removed:
[grains,ebsd.grainId] = calcGrains(ebsd('indexed'),'angle',3*degree);
grains = smooth(grains,5);

% plot the data and the grain boundaries
plot(ebsd('Iron bcc'),ebsd('Iron bcc').orientations,'figSize','large')
hold on
plot(grains.boundary,'linewidth',2)
hold off

%%

job = parentGrainReconstructor(ebsd,grains)

job.p2c = orientation.KurdjumovSachs(job.csParent, job.csChild)

%%


close all
histogram(job.calcGBFit./degree,'BinMethod','sqrt')
xlabel('disorientation angle')

job.calcParent2Child

%%

hold on
histogram(job.calcGBFit./degree,'BinMethod','sqrt')
hold off

%%

% compute the misfit for all child to child grain neighbors
[fit,c2cPairs] = job.calcGBFit;

% select grain boundary segments by grain ids
[gB,pairId] = job.grains.boundary.selectByGrainId(c2cPairs);

% plot the child phase
plot(ebsd('Iron bcc'),ebsd('Iron bcc').orientations,'figSize','large','faceAlpha',0.5)

% and on top of it the boundaries colorized by the misfit
hold on;
% scale fit between 0 and 1 - required for edgeAlpha
plot(gB, 'edgeAlpha', (fit(pairId) ./ degree - 2.5)./2 ,'linewidth',2);
hold off


%%

job.calcVariantGraph

job.clusterGraph

job.votes

%%

plot(job.grains,job.votes.prob(:,1)-job.votes.prob(:,2))

%%

job.calcParentFromVote
%job.calcParentFromVote('minProb',0.5)

%%

plot(job.parentGrains,job.parentGrains.meanOrientation)

%%

% currently not working !
job.selectInteractive

%%

job.mergeSimilar

plot(job.parentGrains,job.parentGrains.meanOrientation)

%%

% currently not working
job.calcGBVotes('reconsiderAll')

%%

plot(job.grains,job.votes.prob(:,1))

%%

job.calcParentFromVote

%%


job.mergeSimilar

job.mergeInclusions

%% calc variant ids

job.calcVariants

job.grains.prop
job.transformedGrains.prop

%%

% associate to each packet id a color and plot
color = ind2color(job.transformedGrains.packetId);
plot(job.transformedGrains,color,'faceAlpha',0.5)

hold on
parentGrains = smooth(job.parentGrains,10);
plot(parentGrains.boundary,'linewidth',3)
hold off
