mtexdata dubna
odf = calcODF(pf)

%%
ori = odf.discreteSample(20000)

h = Miller({1,0,-1,1},{0,2,-2,1},odf.CS);

plotPDF(ori,h,'MarkerSize',1.5,'all','layout',[2 1])

% saveFigure('../pic/scatterV.png')

%%

plotPDF(ori,h,'MarkerSize',1.5,'all','MarkerAlpha',0.1,'layout',[2 1])
%saveFigure('../pic/alphaV.png')

%%

plotPDF(ori,h,'MarkerSize',1.5,'contour','layout',[2 1])
% saveFigure('../pic/contourV.png')

%%

plotPDF(ori,h,'MarkerSize',1.5,'pcolor','layout',[2 1])

mtexColorMap LaboTeX

% saveFigure('../pic/pcolorV.png')


%%

%plot(ori * Miller(1,0,-1,1,odf.CS),'upper','MarkerAlpha',0.1,'MarkerSize',3)

% saveFigure('../pic/alphaV.png')

%%

plot(ori * Miller(1,0,-1,1,odf.CS),'upper','contour','linewidth',2)

% saveFigure('../pic/contourV.png')

%%

plot(ori * Miller(1,0,-1,1,odf.CS),'upper','pcolor')

%%

odf = calcDensity(ori)

%%

h = Miller(1,0,0,ori.CS);
pdf = calcDensity(ori.symmetrise * h)

%%

r = vector3d.Z;
ipdf = calcDensity(inv(ori.symmetrise) * r)

%%

plot(ipdf,'figSize','small')

%saveFigure('../pic/ipdf.png')

%%

contour(ipdf,'figSize','small')
%saveFigure('../pic/ipdfContourf.png')

%%

h = Miller(-1,0,1,1,ori.CS);
pdf = calcDensity(ori.symmetrise * h)

plot3d(pdf)
hold on
arrow3d(xvector)
arrow3d(yvector)
arrow3d(zvector)
hold off
light
%saveFigure('../pic/pdf3d.png')

%%

ori = odf.discreteSample(100)

odf = calcDensity(ori,'halfwidth',15*degree)

plotPDF(odf,h,'layout',[2 1])

%saveFigure('../pic/pdf15.png')
%%

psi = calcKernel(ori)

odf = calcDensity(ori,'kernel',psi)

plotPDF(odf,h,'layout',[2 1])
%saveFigure('../pic/pdfPsi.png')

%%

cs = crystalSymmetry('432')
ss = specimenSymmetry('222')

ori1 = orientation.brass(cs,ss)
ori2 = orientation.goss(cs,ss)

odf = 0.2 * uniformODF(cs,ss) + ...
 0.4 * unimodalODF(ori1,'halfwidth',10*degree) + ...
 0.4 * unimodalODF(ori2,'halfwidth',20*degree);

plot(odf,'sections',9,'layout',[3,3])
annotate(ori1,'displayName','brass')
annotate(ori2)

%saveFigure('../pic/odfSim1.png')

%%

f = fibre.beta(cs,ss)
odf = fibreODF(f,'halfwidth',5*degree)

plot(odf,'sections',9,'layout',[3,3],'resolution',1*degree)

%saveFigure('../pic/odfSim2.png')

%%

plot(odf,'3d')

%saveFigure('../pic/odfSim2_3d.png')


%%


plot(odf,'phi1',(0:10:80)*degree,'layout',[3,3],'resolution',1*degree)

%%

cs = crystalSymmetry('622');

ori1 = orientation.byEuler(30*degree,5*degree,0,cs)
ori2 = orientation.byEuler(30*degree,85*degree,0,cs)

odf = unimodalODF(ori1,'halfwidth',10*degree) + unimodalODF(ori2,'halfwidth',15*degree)

plotSection(odf,'phi2',[0:15:45]*degree,'layout',[4,1])

%saveFigure('../pic/odfSim3_Euler.png')

%%

sigma = linspace(0,60,7)* degree; sigma = sigma(1:end-1);

plotSection(odf,'sigma',sigma,'layout',[3,2])

%saveFigure('../pic/odfSim3_sigma.png')

%%

h = Miller({0,0,0,1},{2,-1,-1,0},{1,-1,0,0},{2,-1,-1,1},cs)
plotPDF(odf,h,'layout',[2 2])

%saveFigure('../pic/odfSim3_pf.png')

%%

pdf = calcPDF(odf,h(1))

%%

r = [vector3d.X,vector3d.Y,vector3d.Z,vector3d(1,1,1)]
plotIPDF(odf,r,'layout',[2,2])

%saveFigure('../pic/odfSim3_ipf.png')

%%

[v, data] = vector3d.load('smiley.csv','columnNames',{'x','y','z','value'})

plot(v, data.value,'layout',[2,1])

%saveFigure('../pic/vecData.png')

%%

sF = interp(v, data.value, 'linear');
plot(sF,'layout',[2,1])

%saveFigure('../pic/vecDataLin.png')

%%

sF = interp(v, data.value, 'harmonic');
plot(sF,'layout',[2,1])

%saveFigure('../pic/vecDataHarm.png')





