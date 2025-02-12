mtexdata dubna

%%

plotx2east

plot(pf{1:2},'layout',[2,1],'MarkerSize',4,'MarkerEdgeColor','gray')
mtexColorMap WhiteJet
mtexColorbar

%saveFigure('../pic/pfRaw.png')

%%


odf = calcODF(pf)

%%

sigma = linspace(0,120*degree,7);
sigma(end) = [];

plot(odf,'sigma',sigma,'layout',[3,2])

% saveFigure('../pic/odfRec.png')


%%

 h = Miller({0,0,0,1},{1,0,-1,0},odf.CS)
plotPDF(odf,h,'layout',[2 1])
%saveFigure('../pic/pfRec.pdf','crop')

%%

pdf = calcPDF(odf,h(2))

%%

[m,pos] = max(pdf,'numLocal',3)


annotate(pos,'parent',gca)
%saveFigure('../pic/pfRecMarked.pdf','crop')

%%


pdf = calcPDF(odf,h(1))

%%

[m,pos] = max(pdf,'numLocal',1)

%%

f = fibre(h(1),pos)

plot(odf,f,'linewidth',2)

xlabel('')
ylabel('')

%saveFigure('../pic/fibre.pdf','crop')

%%

plotx2east

plot(pf{1:2},'layout',[2,1],'MarkerSize',4,'MarkerEdgeColor','gray')
mtexColorMap WhiteJet
mtexColorbar

%saveFigure('../pic/pfRaw2.png')

%%

rot = rotation.rand;

plot(rot*pf{1:2},'layout',[2,1],'MarkerSize',4,'MarkerEdgeColor','gray')
mtexColorMap WhiteJet
mtexColorbar
%saveFigure('../pic/pfRot.png')

%%

plot([1,0.1] * pf{1:2},'layout',[2,1],'MarkerSize',4,'MarkerEdgeColor','gray')
mtexColorMap WhiteJet
mtexColorbar

%saveFigure('../pic/pfScaled.png')

%%

plot([1,0.1] * pf{1:2},'layout',[2,1],'MarkerSize',4,'MarkerEdgeColor','gray')
mtexColorMap WhiteJet
mtexColorbar

%saveFigure('../pic/pfScaled.png')

%%


plot(pf{1:2}.normalize,'layout',[2,1],'MarkerSize',4,'MarkerEdgeColor','gray')
mtexColorMap WhiteJet
mtexColorbar

%saveFigure('../pic/pfNorm.png')

%%

mtexdata geesthacht

% plot imported polefigure
plot(pf{1:2},'layout',[2,1],'MarkerSize',5,'MarkerEdgeColor','gray')
mtexColorMap WhiteJet
mtexColorbar


%saveFigure('../pic/pfGeest.png','crop')

%%

% plot imported polefigure
plot([pf{1},correct(pf{1},'background',pf{2})],'layout',[2,1],'MarkerSize',5,'MarkerEdgeColor','gray')
mtexColorMap WhiteJet
mtexColorbar


%saveFigure('../pic/pfGeestCor.pdf','crop')



%% Splitting and Reordering of Pole Figures
% As we can see the first and the third pole figure complete pole figures
% and the second and the fourth pole figures contain some values for
% background correction. Let us, therefore, split the pole figures into
% these two groups. 

pf_complete = pf({1,3})
pf_background= pf({2,4})

%%
% Actually, it is possible to work with pole figures as with simple numbers.
% E.g. it is possible to add / subtract pole figures. A superposition of
% the first and the third pole figures can be written as

2*pf({1}) + 3*pf({3})


%% Correct pole figure data
%
% In order to correct pole figures for background radiation and
% defocusing one can use the command 
% <PoleFigure.correct.html correct>. In our case the syntax is

pf = correct(pf_complete,'background',pf_background);
plot(pf)

%% Normalize pole figures
%
% Sometimes people want to have normalized pole figures. In the case of
% complete pole figures, this can be simply archived using the command
% <PoleFigure.normalize.html normalize> 

pf_normalized = normalize(pf);
plot(pf_normalized)

%%
% However, in the case of incomplete pole figures, it is well known, that
% the normalization can only by computed from an ODF. Therefore, one has to
% proceed as follows:

% compute an ODF from the pole figure data
odf = calcODF(pf);

% and use it for normalization
pf_normalized = normalize(pf,odf);

plot(pf_normalized)


%% Modify certain pole figure values
% 
% As pole figures are usually experimental data they may contain outliers. In
% order to remove outliers from pole figure data one can use the function
% <PoleFigure.isOutlier.html isOutlier>. Here a simple example:

% Let us add 100 random outliers to the pole figure data
% First we select 100 random positions within the pole figures
ind = randperm(pf.length,100);

% Next we multiply the intensity at these positions by a random value
% between 3 and 4
factor = 3+rand(100,1);
pf(ind).intensities = pf(ind).intensities(:) .* factor;

% Let's check the result
plot(pf)

%%
% check for outliers
condition = pf.isOutlier;

% remove outliers
pf(condition) = [];

% plot the corrected pole figures
plot(pf)

%%
% Sometimes applying the above correction is not suffcient. Then it can
% help to repeat the outlier detection ones again

pf(pf.isOutlier) = [];
plot(pf)

%% Remove certain measurements from the data
% In the same way, as we removed the outlier one can manipulate and delete
% pole figure data by any criteria. Lets, e.g. cap all values that are
% larger than 500. 

% find those values
condition = pf.intensities > 500;

% cap the values in the pole figures
pf(condition).intensities = 500;

plot(pf)


%% Rotate pole figures
% Sometimes it is necessary to rotate the pole figures. In order to do this
% with MTEX one has first to define a rotation, e.e. by

% This defines a rotation around the x-axis about 100 degree
rot = rotation.byAxisAngle(xvector,100*degree);

%%
% Second, the command <PoleFigure_rotate rotate> can be used to rotate the
% pole figure data.
pf_rotated = rotate(pf,rot);
plot(pf_rotated,'antipodal')


%%

plotx2east

plot(pf{1:6},'layout',[3,2],'MarkerSize',4,'MarkerEdgeColor','gray')
mtexColorMap WhiteJet


%saveFigure('../pic/pfRaw6.png')
%%

plotPDF(odf,pf.h(1:6),'upper','layout',[3,2])

%saveFigure('../pic/pfRec6.png')

%%


odf = SantaFe;

plot(odf,'sections',12,'layout',[4 3],'resolution',0.5*degree,...
  'coordinates','off','labels','off','minmax','on','innerPlotSpacing',0)

% saveFigure('../pic/santaFe.png')

%%

h = Miller({1,0,0},{1,1,0},{1,1,1},{2,1,0},{2,1,1},{2,2,1},odf.CS)

pf = calcPoleFigure(odf,h)
figure(2)
plot(pf,'layout',[3,2],'MarkerSize',4,'MarkerEdgeColor','gray')
mtexColorMap WhiteJet
%mtexColorbar

% saveFigure('../pic/santaFePf6.png')

%%

rec2 = calcODF(pf{1:2})
%%
figure(3)
plot(odf,'sections',12,'layout',[4 3],'resolution',1.5*degree,...
  'coordinates','off','labels','off','minmax','on','innerPlotSpacing',0)

% saveFigure('../pic/santaFeRec2.png')

%%

rec = calcODF(pf,'noGhostCorrection')
%%
figure(3)
plot(odf,'sections',12,'layout',[4 3],'resolution',1.5*degree,...
  'coordinates','off','labels','off','minmax','on','innerPlotSpacing',0)

% saveFigure('../pic/santaFeRecNoGhost.png')

%%

plotFourier(SantaFe,'DisplayName','odf',...
  'Marker','s','MarkerSize',8,'MarkerFaceColor','k','linestyle','none','MarkerEdgeColor','k')

hold on
plotFourier(rec,'linewidth',3,'DisplayName','rec','linestyle','none','MarkerSize',12)
plotFourier(rec2,'linewidth',2,'DisplayName','rec2','linestyle','none','MarkerSize',14)
hold off
legend show
xlim([0,17])

set(gcf,'Renderer','painters')

% saveFigure('../pic/santaFeSpectra.pdf','crop')
