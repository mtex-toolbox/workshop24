

%%%

cS = crystalShape.quartz;
mori = orientation.byAxisAngle(Miller(1,0,-1,1,cS.CS),87*degree);

plot(cS,'colored')
hold on
plot(mori*cS, 'colored','FaceAlpha',0.5)
hold off
legend off

%saveFigure('../pic/japanTwin.png')

%%

mori.angle('max') /degree

%%

v = mori.axis('max')


%% japan twins at xi -negative trigonal bipyramid

mori = orientation.byAxisAngle(Miller({1 1 -2 2},cS.CS),180*degree); 

plot(cS,'colored')
hold on
plot(mori*cS, 'colored','FaceAlpha',0.5)
hold off
legend off

%%

mori = orientation.byAxisAngle(Miller({1 1 -2 2},cS.CS),180*degree); 
angle(cS.CS.cAxis, mori * cS.CS.cAxis)./degree


%%

mtexdata alphaBetaTitanium

beta2alpha = orientation.Burger(ebsd('Ti (beta)').CS,ebsd('Ti (alpha)').CS);

round2Miller(beta2alpha)

%%
csBeta = ebsd('Ti (beta)').CS;
csAlpha = ebsd('Ti (alpha)').CS;
oBeta = orientation.byEuler(0,0,0,csBeta)
r = oBeta * Miller(1,1,0,csBeta)

%%


oAlpha = oBeta * inv(beta2alpha)

oAlpha.symmetrise

%%

oAlpha = orientation.byEuler(135*degree,90*degree,55*degree,csAlpha)
inv(oAlpha) * r

%%

b2a = inv(oAlpha) * oBeta

%%

b2a * Miller(1,1,0,csBeta)


inv(b2a) * Miller(0,0,0,1,csAlpha)

%%

round2Miller(b2a)

%%

plotIPDF(b2a.symmetrise,Miller(0,0,0,1,csAlpha),'MarkerSize',10)

%%

oriBeta = orientation.id(csBeta)

oriAlpha = variants(beta2alpha,oriBeta)

hBeta = Miller({1,1,0},{1,-1,1},csBeta);
hAlpha1 = Miller(0,0,0,1,csAlpha);
hAlpha2 = Miller(-2,1,1,0,csAlpha,'uvw')
hBeta2 = Miller({1,-1,1},csBeta,'uvw');

for k = 1:12

  plotPDF(oriBeta,hBeta,'figSize','large','layout',[2 1])
  
  nextAxis(1)
  %xlabel(char(hAlpha1,'latex'),'Interpreter','latex')
  mtexTitle([char(hBeta(1),'latex') '  ' char(hAlpha1,'latex')])
  hold on
  plot(oriAlpha(k) * symmetrise(hAlpha1),...
    'MarkerFaceColor','none','MarkerSize',15,'LineWidth',4,...
    'MarkerEdgeColor',ind2color(1+k))
  nextAxis(2)
  hold on
  %xlabel(char(hAlpha2,'latex'),'Interpreter','latex')
  mtexTitle([char(hBeta2,'latex') '  ' char(hAlpha2,'latex')])
  plot(oriAlpha(k) * hAlpha2.symmetrise,...
    'MarkerFaceColor','none','MarkerSize',15,'LineWidth',4,...
    'MarkerEdgeColor',ind2color(1+k))
  hold off
  set(gcf,'Renderer','painters')
  drawNow(gcm)
  %saveFigure(['../pic/varPf' int2str(k) '.pdf'],'crop')

end

%%

plotPDF(oriBeta,hBeta,'figSize','large','layout',[2 1])
  
nextAxis(1)
mtexTitle([char(hBeta(1),'latex') '  ' char(hAlpha1,'latex')])

hold on
  plot(oriAlpha * symmetrise(hAlpha1),...
    'MarkerFaceColor','none','MarkerSize',15,'LineWidth',3,...
    'MarkerEdgeColor','red')
  
  nextAxis(2)
  hold on
  mtexTitle([char(hBeta2,'latex') '  ' char(hAlpha2,'latex')])
  
  plot(oriAlpha * hAlpha2.symmetrise,...
    'MarkerFaceColor','none','MarkerSize',15,'LineWidth',3,...
    'MarkerEdgeColor','red')
  hold off
  set(gcf,'Renderer','painters')
  drawNow(gcm)
  
  saveFigure(['../pic/varPf.pdf'],'crop')

%%

calcVariantId(oriBeta,oriAlpha,beta2alpha)
  
%%

oriA = orientation.rand(1000,3,csAlpha);

[~,fit] = calcParent(oriA,beta2alpha);

hold on
histogram(fit./degree)
xlabel('fit in degree')

%saveFigure('../pic/fitHist3.pdf','crop')

%%

n = 4;

oriB = orientation.rand(1000,csBeta);

b2aV = variants(beta2alpha);

vId = randi(12,length(oriB),n);

oriA = rotation.rand(1000,n,'maxAngle',2*degree) .* ...
  repmat(oriB,1,n) .* inv(b2aV(vId));

[~,fit] = calcParent(oriA,beta2alpha,'numFit',2);

histogram(fit(:,1)./degree)
hold on
histogram(fit(:,2)./degree)
hold off
xlabel('fit in degree')

%saveFigure('../pic/fit3.pdf','crop')




%%

mtexdata alu
ebsd = ebsd.gridify

ori = ebsd.orientations;

%plot(ebsd,1e-1 + norm(ebsd.gradientX)./degree)
plot(ebsd,1e-1 + norm(ebsd.gradientY)./degree)

setColorRange([0,12])
%mtexColorbar
set(gca,'ColorScale','log')

%saveFigure('../pic/gradientX.png')
%saveFigure('../pic/gradientY.png')
%%

[grains, ebsd.grainId] = calcGrains(ebsd('indexed'))
grains =smooth(grains,5)
%%
plot(ebsd,ebsd.orientations)

hold on
plot(grains.boundary('indexed'),'linewidth',7)
plot(grains.boundary,grains.boundary.misorientation.angle./degree,'linewidth',5)
hold off
mtexColorbar

set(gcf,'Renderer','opengl')

%saveFigure('../pic/boundaryAngle.png')

%%

oriRef = grains('id',ebsd('indexed').grainId).meanOrientation

plot(ebsd('indexed'),angle(ebsd('indexed').orientations,oriRef)./degree)
mtexColorMap LaboTeX
hold on
plot(grains.boundary,'lineWidth',2)
hold off
mtexColorbar

%saveFigure('../pic/mis2mean.png')
%%

plot(grains,grains.GOS./degree,'lineWidth',2)
mtexColorMap LaboTeX
mtexColorbar

%saveFigure('../pic/gos.png')

%%

mtexdata twins

[grains, ebsd.grainId] = calcGrains(ebsd('indexed'))

%%

histogram(grains.boundary.misorientation.angle./degree)
%saveFigure('../pic/angleHist.pdf','crop')

%%
plotAngleDistribution(grains.boundary.misorientation)

%%

  bmdf = calcDensity(mori,'halfwidth',5*degree)
  plotAngleDistribution(bmdf)

%%

odf = calcODF(mori)
mdf = calcMDF(odf)

%%

plotAngleDistribution(grains.boundary.misorientation)
xlim([0,oR.maxAngle ./degree])

hold on
plotAngleDistribution(odf.CS,'displayName','uniform')
hold off

legend show

%saveFigure('../pic/angleUMDF.pdf','crop')



%%

hold on
plotAngleDistribution(mdf,'displayName','umdf')
hold off
legend show

%saveFigure('../pic/angleMDF.pdf','crop')

%%

hold on
plotAngleDistribution(odf.CS,'displayName','umdf')
hold off

%%
mori = grains.boundary.misorientation;
bmdf = calcDensity(mori,'halfwidth',5*degree)
hold on
plotAngleDistribution(bmdf)
hold off

%saveFigure('../pic/angleFull.pdf','crop')

%%  
  
oR =  odf.CS.fundamentalRegion

%%

plotAxisDistribution(mori)

%saveFigure('../pic/axisScatter.png')

%%

plotAxisDistribution(mori,'contourf')
mtexColorbar
mtexColorMap LaboTeX
%saveFigure('../pic/axisSmooth.png')

%%
plotx2south
plotAxisDistribution(uniformODF(mori.CS,mori.CS),'antipodal')
mtexColorMap LaboTeX
%saveFigure('../pic/axisCS.png')

%%

%plotx2south
plotAxisDistribution(mdf,mori.CS.fundamentalSector)
mtexColorMap LaboTeX
%saveFigure('../pic/axisMDF.png')

%%

plotAxisDistribution(mori.CS,mori.CS.fundamentalSector)

mtexColorMap LaboTeX
hold on
plotAxisDistribution(mori,'MarkerFaceColor','black','MarkerSize',3)
hold off
