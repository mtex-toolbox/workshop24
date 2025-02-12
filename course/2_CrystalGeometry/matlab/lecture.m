
abc = [5.2, 5.2, 5.3 ]
abg = [90, 99.25 90] * degree;
cs = crystalSymmetry ('1' , abc , abg )

%%

d = Miller(1, 1, 0 ,cs ,   'uvw')




%%

cs = crystalSymmetry.load('Mg-Magnesium.cif')

%cs.axes(:,3) = 5*cs.axes(:,3)

cS = crystalShape.hex(cs)

plot(cS,'colored')

%%


cs = crystalSymmetry.load('quartz.cif')

cS = crystalShape.quartz;

plot(cS,'colored')
legend off

%saveFigure('../pic/quartzMillerColor.png')



%%
c = Miller({0,0,0,1},cs);  % c-axis
m = Miller({1,0,-1,0},cs);  % hexagonal prism
r = Miller({1,0,-1,1},cs);  % positive rhomboedron, usally bigger then z
z = Miller({0,1,-1,1},cs);  % negative rhomboedron
s1 = Miller({2,-1,-1,1},cs);% left tridiagonal bipyramid
s2 = Miller({1,1,-2,1},cs); % right tridiagonal bipyramid
x1 = Miller({6,-1,-5,1},cs);% left positive Trapezohedron
x2 = Miller({5,1,-6,1},cs); % right positive Trapezohedron

h = [m,r,z,s2,x2,c];



sR = sphericalRegion('maxTheta',pi/2,cs.plotOptions{:})
proj = 'eangle';
%proj = 'edist';
%proj = 'orthographic';
%proj = 'earea'; 
%proj = 'gnonomic'; sR = sphericalRegion('maxTheta',s1.theta,cs.plotOptions{:})



newSphericalPlot(sR,'projection',proj,'grid','on','figSize','medium','grid_Res',10*degree)
hold on
for i = 1:length(h)
  plot(h(i).symmetrise,'labeled','upper','backgroundColor','White')
end
hold off

%saveFigure(['../pic/' proj '.pdf']);

%%

clf
for i = 1:length(h)
  plot(h(i).symmetrise,'3d',...
    'displayName',char(h(i)))
  hold on
end
hold off
%legend show Location eastoutside
%saveFigure(['../pic/proj3d.png']);

%%




