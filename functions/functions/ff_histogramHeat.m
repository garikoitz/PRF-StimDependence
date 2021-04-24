function ff_histogramHeat(x, y, maxValueX, maxValueY, numHistBins,cmapValuesHist,fov,roiName)
% ff_histogramHeat(x, y, maxValueX, maxValueY, numHistBins) 
%
% Makes a heat map!
% INPUTS
% x:            the vector of x values. assuming row vec!
% y:            the vector of y values. assuming row vec!
% maxValue:     limits of the x and y axis
% numHistBins   something around 50

% GLU: used numHistBins as the radius of scatplot now, it was unused

%% 
minValueX = 0; %-2.3; % 0
minValueY = 0; %-2.3; % 0

axisLimsX = [minValueX maxValueX]; 
axisLimsY = [minValueY maxValueY]; 
% Ctrs{1} = linspace(minValueX, maxValueX, numHistBins);
% Ctrs{2} = linspace(minValueY, maxValueY, numHistBins);

% the histogram!
% figure
% hist3([x' y'],'Ctrs', Ctrs)
% set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')
% hold on


% New GLU: 2D plot, simpler and more control
% scatplot(x,y,method,radius,N,n,po,ms,colormap)
scatplot(x,y,'squares',numHistBins,100,5,1,10,cmapValuesHist);
set(gca,'color','k');
set(gca, 'xlim', axisLimsX);
set(gca, 'ylim', axisLimsY);
set(gca, 'FontSize', 14);
axis square;
hold on
% identityLine goes above everything else so that it can be seen
% iColor = [0 1 1]; % cyan
% iColor = [1 1 0]; % yellow
iColor = [.2 .5 .2]; % blue
npoints = 100; 
% maxZ = max(get(gca, 'ZLim'));
% zVec = maxZ*ones(1, npoints); 

% plot3(linspace(minValueX, maxValueX, npoints), linspace(minValueY, maxValueY, npoints), zVec, ...
%     '--', 'Color', iColor, 'LineWidth',2)
% plot3(linspace(minValueX+fov, maxValueX, npoints), linspace(minValueY, maxValueY-fov, npoints), zVec, ...
%     '--', 'Color', iColor, 'LineWidth',1.5)
% plot3(linspace(minValueX, maxValueX-fov, npoints), linspace(minValueY+fov, maxValueY, npoints), zVec, ...
%     '--', 'Color', iColor, 'LineWidth',1.5)
% jbfill([0.8, length(upperVals)+.2],[-2,-2],[2,2],[.1,.1,.1],[.2,.2,.2],0,0.1);



line([minValueX, maxValueX], [minValueY, maxValueY], ...
      'LineStyle','--', 'Color', iColor, 'LineWidth',2)
% hline(linspace(minValueX, maxValueX-fov, npoints), ...
%       linspace(minValueY+fov, maxValueY, npoints), zVec, ...
%     '--', 'Color', iColor, 'LineWidth',1.5)
jbfill([minValueX,minValueX+fov, maxValueX-fov, maxValueX], ...
       [minValueY+fov,minValueY+fov+fov,maxValueY,maxValueY],...
       [minValueY,minValueY,maxValueY-fov-fov,maxValueY-fov], ...
       iColor,iColor,0,0.4);
line([minValueX, maxValueX], [5,5], ...
      'LineStyle','-', 'Color', 'r', 'LineWidth',1)
roiName = [strrep(ff_stringRemove(roiName, 'WangAtlas_'),'_','\_')];
text(maxValueX-(maxValueX-minValueX)*(1/4), minValueY+fov,roiName,...
     'Color','w','FontWeight','Bold','FontSize',14)
colormap(cmapValuesHist); 
c = colorbar;
set(c, 'location', 'eastoutside')








% color things         
% matlab has funky behavior where the size of this influences the size of 
% all future colorbars...
% colormap(zeros(64,3)); 
% cmapValuesHist = colormap('pink');


% caxis([0 maxZ]);
% caxis([0 maxZ/5])
        
end