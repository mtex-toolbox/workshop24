
figure(1)
ori = orientation.rand(10000,2,fcc2bcc.CS);

[~,fit] = calcParent(ori,fcc2bcc,'numFit',2);

fit = fit(fit(:,1)<10*degree,:);


histogram(fit(:,1) ./ degree,100,'Normalization','probability')

hold on
histogram(fit(:,2) ./ degree,100,'Normalization','probability')
hold off
xlabel('fit in degree')

xlim([0,10])

%%

figure(2)
scatter(fit(:,1),fit(:,2),'MarkerEdgeAlpha',0.05)



%%

figure(2)
ori = orientation.rand(10000,3,fcc2bcc.CS);

[~,fit] = calcParent(ori,fcc2bcc,'numFit',2);

%fit = fit(fit(:,1)<10*degree,:);

histogram(fit(:,1) ./ degree,100)

hold on
histogram(fit(:,2) ./ degree,100)
hold off
xlabel('fit in degree')
xlim([0,9])
ylim([0,50])

%%

[~, fit] = calcParent(grains(grainPairs).meanOrientation,fcc2bcc,'numFit',2);

%%

histogram(fit(:,1) ./ degree,200)
hold on
histogram(fit(:,2) ./ degree,200)
hold off
xlabel('fit in degree')
xlim([0,5])
set(gca,'YTickLabel',{})

%saveFigure('../pic/histGB.pdf')

%%


tP = grains.triplePoints;
tPori = grains(tP.grainId).meanOrientation;

% compute the misfit to a common parent orientation
[~, fit2] = calcParent(tPori,fcc2bcc,'id','threshold',10*degree,'numFit',2);
%%

histogram(fit2(:,1) ./ degree,200)
hold on
histogram(fit2(:,2) ./ degree,200)
hold off
xlabel('fit in degree')
xlim([0,5])

%saveFigure('../pic/histTP.pdf')




