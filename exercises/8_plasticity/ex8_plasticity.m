%%                        MTEX Workshop 2024
%
%%                   Exercise 9 - Plasticity


%% Magnesium

cs = crystalSymmetry.load('Mg-Magnesium.cif')
cs = cs.properGroup;

% the initial ODF
odf = fibreODF(cs.cAxis, vector3d.Z,'halfwidth',20*degree);
odf = uniformODF(cs)
odf = SO3FunHarmonic(odf)


sS = [slipSystem.basal(cs,1),...
  slipSystem.prismatic2A(cs,15),...
  slipSystem.pyramidalCA(cs,10),...
  slipSystem.twinC1(cs,100)];

% consider all symmetrically equivalent slip systems
sS = sS.symmetrise;

eps = 0.35 * strainTensor(diag([1 -0.5 -0.5]));

%% the Taylor factor

[M,~,W] = calcTaylor(eps, sS)

% for some reason we need some smoothing of the vector field
psi = SO3DeLaValleePoussinKernel('halfwidth',5*degree);
W.SO3F = W.SO3F.conv(psi);

%% the orientation dependence of the Taylor factor

plot(M,'sigma',[0 30]*degree)

%%

plotPDF(M,Miller(0,0,0,1,cs))
mtexColorbar

%%

plotIPDF(M,[xvector,yvector,zvector])

%% the spin tensor

% its divergence
plot(div(W),'sigma',[0,30]*degree,'resolution',0.5*degree,'noGrid')
mtexColorMap blue2red

% the spin tensor as an vector field
hold on
plot(W,'add2all','noGrid','color','black','resolution',5*degree)
setColorRange([-3,3])
hold off


%% texture evolution - the traditional route

ori = odf.discreteSample(10000);

ori1 = doEulerStep(@(ori) -calcSpin(ori,eps,sS),ori,5)
odf1 = calcDensity(ori1,'halfwidth',2.5*degree);

figure(1)
plot(odf1,'sigma')
mtexColorbar

%% texture evolution - a bit quicker

ori2 = doEulerStep(-W,ori,5)
odf2 = calcDensity(ori2,'halfwidth',2.5*degree);

figure(2)
plot(odf2,'sigma')
mtexColorbar

%% texture evolution - without ODF discretisation

odf3 = doEulerStep(-W,odf,5)

figure(3)
plot(odf3,'sigma')
mtexColorbar


%% FCC material

cs = crystalSymmetry('432');

% define a family of slip systems
sS = slipSystem.fcc(cs);
sS = sS.symmetrise;

% 30 percent plane strain
q = 0;
eps = 0.3 * strainTensor(diag([1 -q -(1-q)]));

odf = uniformODF(cs);

%% 
[~,~,W] = calcTaylor(eps, sS)
psi = SO3DeLaValleePoussinKernel('halfwidth',5*degree);
W.SO3F = W.SO3F.conv(psi);

%% the traditional route

% define some random orientations
ori = odf.discreteSample(1e5);

ori1 = doEulerStep(@(ori) calcSpin(ori,eps,sS),ori,5)
odf1 = calcDensity(ori1,'halfwidth',2.5*degree);

figure(1)
plot(odf1,'sigma')


%% route 2

ori2 = doEulerStep(-W,ori,5)

odf2 = calcDensity(ori2,'halfwidth',5*degree);

figure(2)
plot(odf2,'sigma')

%% route 3

odf3 = doEulerStep(-W,odf,10)

figure(3)
plot(odf3,'sigma')
mtexColorbar


%% old route 

ori4 = ori;

numIter = 10;
progress(0,numIter);
for sas=1:numIter

  % compute the Taylor factors and the orientation gradients
  WLocal = spinTensor(W.eval(ori4).').';

  % rotate the individual orientations
  ori4 = ori4 .* orientation(-WLocal./numIter);
  progress(sas,numIter);

end

odf4 = calcDensity(ori4,'halfwidth',5*degree);

figure(4)
plot(odf4,'sigma')




%% Single Slip Model

csOli = crystalSymmetry('222',[4.779 10.277 5.995],'mineral','olivine');

sSOli = slipSystem(Miller({1,0,0},{1,0,0},{0,0,1},csOli,'uvw'),...
  Miller({0,1,0},{0,0,1},{0,1,0},csOli,'hkl'))

S = sSOli.deformationTensor;

E = 0.3 * strainRateTensor([1 0 0; 0 0 0; 0 0 -1]);

Omega = @(ori) -SO3TangentVector(spinTensor(((ori * S(2)) : E) .* (ori * S(2))));
% ori = orientation.rand(csOli)
% 0.5 * EinsteinSum(tensor.leviCivita,[1 -1 -2],(ori * S(2) : E) .* (S(2)),[-1 -2])

Omega = SO3VectorFieldHarmonic.quadrature(Omega,csOli);


%% route - exact solution

figure(1)
odf1 = SO3FunSBF(sSOli(2),E)
plotSection(odf1,'sigma')
mtexColorbar

%% route 2 - Euler equation

% the starting ODF
odf0 = uniformODF(csOli);

odf2 = doEulerStep(2*Omega,odf0,40)

figure(2)
plot(odf2,'sigma')
mtexColorbar


%% route 3

numIter = 1;
res = 2.5*degree;

%ori0 = equispacedSO3Grid(csOli,'resolution',res);

ori0 = odf0.discreteSample(10000);

ori = doEulerStep(2*Omega,ori0,numIter);

odf3 = calcDensity(ori,'halfwidth',res*4);

figure(3)
plot(odf3,'sigma')
mtexColorbar

%%

plotSection(ori0,'sigma',0,'noGrid','all')
%saveFigure('~/presentations/23/manchester/pic/ori0.pdf')

%%

plotSection(ori,'sigma',0,'noGrid','all')
%saveFigure('~/presentations/23/manchester/pic/ori1.pdf')

%%

% We may visualize the orientation depedence of the spin tensor by plotting
% its divergence in sigma sections and on top of it the spin tensors as a
% quiver plot

plotSection(div(Omega),'sigma',0,'noGrid','contourf','facealpha',0.5)
mtexColorMap blue2red
mtexColorbar

hold on
plot(Omega,'add2all','linewidth',2,'color','k','resolution',7.5*degree)
hold off

%saveFigure('~/presentations/23/manchester/pic/Omega.pdf')

%%

plotSection(odf1,'sigma',0,'noGrid','contourf')
mtexColorbar
%saveFigure('~/presentations/23/manchester/pic/odfTrue.pdf')

%%

numIter = 5;
for k = 1:numIter
  odf2 = doEulerStep(2*Omega*k/numIter,odf0,k);
  plotSection(odf2,'sigma',0,'noGrid','contourf')
  mtexColorbar
  saveFigure(['~/presentations/23/manchester/pic/odfE' num2str(k) '.pdf'])
end

%% 




%%

function W = calcSpin(ori,eps,sS)
[~,~,W] = calcTaylor( inv(ori) .* eps, sS);
end
