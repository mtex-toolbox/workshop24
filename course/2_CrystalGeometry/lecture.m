cs1 = crystalSymmetry('622')

ori = orientation.rand(cs)




%%

cs2 = crystalSymmetry('622','x||a')

%%

ori.transformReferenceFrame(cs2)



cs1 == cs2

%%

h = Miller(1,2,3,cs1)

h.transformReferenceFrame(cs2)

%%

sF = S2FunHarmonicSym(1,cs1);

odf = unimodalODF(cs1)

odf.transformReferenceFrame

%%

ori = orientation.byEuler(10*degree,20*degree,30*degree,cs1);


h = Miller(1,1,-2,0,cs1);

ori * h

% 1 1 ->  0.702179  1.77652 0.592396
% 1 2 ->  0.702179  1.77652 0.592396
% 2 1 -> -0.326352   1.85083   0.68404
% 2 2 -> -0.326352   1.85083   0.68404


%%

ori = orientation.byEuler(10*degree,20*degree,15*degree,cs1);

plotSection(ori)


%%

cs = crystalSymmetry('432')


ori = orientation.rand(10,cs)

plot(ori)


%%

plot(ori.symmetrise)

%%

for i=1:length(ori)

  plot(ori(i).symmetrise,'filled','MarkerEdgeColor',ind2color(i))
  hold on

end
hold off

%%

plot(ori,'axisAngle','filled')

%%
for i=1:length(ori)
  
  plotSection(ori(i),'MarkerFaceColor',ind2color(i))
  hold on

end
hold off

%%

for i=1:length(ori)
  
  plotSection(ori(i),'MarkerFaceColor',ind2color(i),'sigma','sections',12)
  hold on

end
hold off

%%

plotSection(ori,ind2color(1:length(ori)),'axisAngle',(5:5:60)*degree)


%%

ori = orientation.rand(100000);

plotSection(ori,'MarkerFaceAlpha',0.5,'all','sections',1)


%%




