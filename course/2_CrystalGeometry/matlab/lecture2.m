
cs  = crystalSymmetry('622')
euler = [10,20,30]*degree;  
ori = orientation.byEuler(euler,cs)
h   = Miller(1,1,-2,0,cs)

%%

inv(ori)

%%


cs = crystalSymmetry('432')
ss = specimenSymmetry('222')
o1 = orientation.goss(cs,ss)
o2 = orientation.brass(cs,ss)

f = fibre(o1,o2)

%%

f = fibre(o1,o2,'full')

%%



cs_x2a = crystalSymmetry('622','X||a*','Z||c');

% visualize the resuls
plot(cs_x2a,'figSize','small')
annotate(cs_x2a.aAxis,'MarkerFaceColor','r','label','a','backgroundColor','w')
annotate(cs_x2a.bAxis,'MarkerFaceColor','r','label','b','backgroundColor','w')
annotate(-vector3d.Y,'MarkerFaceColor','green','label','-y','backgroundColor','w')
annotate(-vector3d.X,'MarkerFaceColor','green','label','-x','backgroundColor','w')

%saveFigure('../pic/622xastar.pdf')


%%
% In contrast the following command alignes the $\vec y$ axes to the $\vec
% a$ axes and the $\vec z$ axes to the $\vec c$ axes.

cs_y2a = crystalSymmetry('622','x||a','Z||c');
plot(cs_y2a,'figSize','small')
annotate(cs_y2a.aAxis,'MarkerFaceColor','r','label','a','backgroundColor','w')
annotate(cs_y2a.bAxis,'MarkerFaceColor','r','label','b','backgroundColor','w')
annotate(-vector3d.Y,'MarkerFaceColor','green','label','-y','backgroundColor','w')
annotate(-vector3d.X,'MarkerFaceColor','green','label','-x','backgroundColor','w')

%saveFigure('../pic/622xa.pdf')

%%
% The only difference between the above two plots is the position of the
% $\vec x$ and $\vec y$ axes. The reason is that visualizations relative to
% the crystal reference system, e.g., inverse pole figures, are in MTEX
% aligned on the screen according to the b-axis.
%
% This on secreen alignment can be easily modified by

% change on screen alignment
plota2east

% redo last plot
plot(cs_y2a,'figSize','small')
annotate(cs_y2a.aAxis,'MarkerFaceColor','r','label','a','backgroundColor','w')
annotate(cs_y2a.bAxis,'MarkerFaceColor','r','label','b','backgroundColor','w')
annotate(-vector3d.Y,'MarkerFaceColor','green','label','-y','backgroundColor','w')
annotate(-vector3d.X,'MarkerFaceColor','green','label','-x','backgroundColor','w')

% set old default back
plotb2east

%%
% It should be stressed that the alignment between the Eucledean crystal
% axes $\vec x$, $\vec y$, $\vec z$ and the crystallographic axes $\vec a$,
% $\vec b$ and $\vec c$ is crucial for many computations. The difference
% between both setups becomes more vsible if we plot crystal shapes in the
% $\vec x$, $\vec y$, $\vec z$ coordinate system

cS_x2a = crystalShape.quartz(cs_x2a);

close all
figure(1)
plot(cS_x2a,'figSize','small','colored')
hold on
arrow3d(0.5*[xvector,yvector,zvector],'labeled')
hold off

%%

cS_y2a = crystalShape.quartz(cs_y2a);

figure(2)
plot(cS_y2a,'figSize','small')
hold on
arrow3d(0.5*[xvector,yvector,zvector],'labeled')
hold off


%%
% Most important is the difference if Euler angles are used to describe
% orientation. Lets consider the following two orientations

ori_x2a = orientation.byEuler(0,0,0,cs_x2a)
ori_y2a = orientation.byEuler(0,0,0,cs_y2a)

%%
% and visualize them in a pole figure. 

newMtexFigure('innerPlotSpacing',20,'figSize','small')
plotPDF(ori_x2a,Miller(1,0,0,cs_x2a),'MarkerSize',20)
annotate([xvector,yvector],'label',{'x','y'},'backgroundColor','w')
nextAxis
plotPDF(ori_y2a,Miller(1,0,0,cs_y2a),'MarkerSize',20)
annotate([xvector,yvector],'label',{'x','y'},'backgroundColor','w')

%%
% We observe that both pole figures are rotated with respect to each other
% by 30 degree. Indeed computing the misorientation angle between both
% orientations gives us

angle(ori_x2a, ori_y2a) ./ degree

%%
% In many cases MTEX automatically recognizes different setups and corrects
% for this. In order to manually transform orientations or tensors from one
% reference frame into another reference frame one might use the command
% <orientation.transformReferenceFrame.html transformReferenceFrame>. The
% following command transfroms the reference frame of orientation |ori_y2a|
% into the reference frame |cs_x2a|

ori_x2a.transformReferenceFrame(cs_y2a)


%% Triclinic and monoclinic symmetries
%
% In triclinic and monoclinic symmetries even more different setups are
% used. As two perpedicular crystal axes are required to align with $\vec
% x$, $\vec y$ or $\vec z$ one ussually chooses one crystal axis from the
% direct coordinate system, i.e., $\vec a$, $\vec b$ or $\vec c$, and the
% second crystal axis from the reciprocal axes $\vec a^*$, $\vec b^*$ or
% $\vec c^*$. Typical examples for such setups are

cs = crystalSymmetry('-1', [8.290 12.966 7.151], [91.18 116.31 90.14]*degree,...
  'x||a*','y||b', 'mineral','An0 Albite 2016')

%%
% or

cs = crystalSymmetry('-1', [8.290 12.966 7.151], [91.18 116.31 90.14]*degree,...
  'x||a','c||c*', 'mineral','An0 Albite 2016')

%%

inv(ori) * vector3d.Z

%%

cs = crystalSymmetry('622')
ori = orientation.byEuler(10*degree,20*degree,30*degree,cs)
h = Miller(1,1,-2,0,cs)
ori * h

%%

inv(ori)

%%

ori = orientation.rand(100,cs);

plot(ori)

%saveFigure('../pic/Euler3d.pdf')

%%

plot(ori,'axisAngle','linewidth',2)
axis off
%saveFigure('../pic/AxisAngle3d.png')

%%

plot(ori,'rodrigues','linewidth',2)
axis off
%saveFigure('../pic/rodrigues3d.png')

%%

plotSection(ori)
%saveFigure('../pic/EulerSection.pdf')

%%

plotSection(ori,'axisAngle',(10:10:80)*degree)
%saveFigure('../pic/AxisAngleSection.pdf')

%%

plotSection(ori,'sigma','sections',12)
%saveFigure('../pic/SigmaSection.pdf')

%%

h = Miller({1,0,-1,0},{1,1,-2,0},{1,0,-1,2},{0,0,0,1},cs)
plotPDF(ori,h,'MarkerSize',5)
%saveFigure('../pic/PFRandom.pdf')

%%

r = [xvector,yvector,zvector]
plotIPDF(ori,r,'noLabel')
annotate(Miller(0,0,0,1,cs),'labeled')
annotate(Miller(-1,2,-1,0,cs),'labeled')
annotate(Miller(-2,1,1,0,cs),'labeled','backgroundColor','w')
%saveFigure('../pic/IPFRandom.pdf')

%% 

plotx2east

odf = SO3Fun.dubna
cs = odf.CS;
h = Miller({0,0,0,1},{1,0,-1,0},{1,0,-1,1},{1,1,-2,1},cs)

plotPDF(odf,h)


%%


tic;[~,ori] = max(odf);toc
%ori = odf.discreteSample(1000)


hold on
plotPDF(ori,h)
hold off


%%


odf.calcComponents

