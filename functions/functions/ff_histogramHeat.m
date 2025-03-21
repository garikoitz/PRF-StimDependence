function c = ff_histogramHeat(x, y, minmaxX, minmaxY, numHistBins,...
                             cmapValuesHist,fov,roiName,fieldName,fontsize, cutoff)
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
minValueX = minmaxX(1);
maxValueX = minmaxX(2);

minValueY = minmaxY(1);
maxValueY = minmaxY(2);


axisLimsX = minmaxX; 
axisLimsY = minmaxY; 
% Ctrs{1} = linspace(minValueX, maxValueX, numHistBins);
% Ctrs{2} = linspace(minValueY, maxValueY, numHistBins);

% the histogram!
% figure
% hist3([x' y'],'Ctrs', Ctrs)
% set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')
% hold on


% New GLU: 2D plot, simpler and more control
% scatplot(x,y,method,radius,N,n,po,ms,colormap)
outdata = scatplot(x,y,'squares',numHistBins,100,5,1,10,cmapValuesHist, cutoff);
set(gca,'color','k');
set(gca, 'xlim', axisLimsX);
set(gca, 'ylim', axisLimsY);
set(gca, 'FontSize', fontsize);
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

% Plot a big red dot to see where the max value of the density is
% plot(outdata.upper_max_dens_coord_x, ...
%      outdata.upper_max_dens_coord_y, 'ro', 'LineWidth', 3)
% plot(outdata.lower_max_dens_coord_x, ...
%      outdata.lower_max_dens_coord_y, 'ro', 'LineWidth', 3)
% plot(outdata.global_max_dens_coord_x, ...
%      outdata.global_max_dens_coord_y, 'r+', 'LineWidth', 3)
% 
% line([minValueX, maxValueX], [minValueY, maxValueY], ...
%       'LineStyle','--', 'Color', iColor, 'LineWidth',2)
% hline(linspace(minValueX, maxValueX-fov, npoints), ...
%       linspace(minValueY+fov, maxValueY, npoints), zVec, ...
%     '--', 'Color', iColor, 'LineWidth',1.5)
jbfill([minValueX,minValueX+fov, maxValueX-fov, maxValueX], ...
       [minValueY+fov,minValueY+fov+fov,maxValueY,maxValueY],...
       [minValueY,minValueY,maxValueY-fov-fov,maxValueY-fov], ...
       iColor,iColor,0,0.4);
% Guess hemi
hemi = 'L ';
if contains(lower(roiName),'right'); hemi = 'R '; end
if strcmp(hemi,'L ')
    roiName = strrep(ff_stringRemove(roiName, 'WangAtlas_'),'_left','');
    roiName = strrep(roiName,'l','');
end
if strcmp(hemi,'R ')
    roiName = strrep(ff_stringRemove(roiName, 'WangAtlas_'),'_right','');
    roiName = strrep(roiName,'r','');
end
text(maxValueX-(maxValueX-minValueX)*(1/2), minValueY+fov/2, ...
    sprintf('%s %s (N:%i)',hemi, roiName, size(x,2)) , ...
     'Color','w','FontWeight','Bold','FontSize',fontsize+4)
colormap(cmapValuesHist); 
c = colorbar;
set(c, 'location', 'eastoutside')



switch fieldName
    case {'ecc'}
        % line([minValueX, maxValueX], [cutoff,cutoff], 'LineStyle','-', 'Color', 'r', 'LineWidth',1)
        % Here we can do the ecc tests in BarData1 and BarData2
        % Calculate standard dev
        S = std(x-y);
        % Cohens d
        fovealInd = (y <= cutoff);
        dfoveal = computeCohen_d(y(fovealInd), x(fovealInd),'paired');
        dperiph = computeCohen_d(y(~fovealInd), x(~fovealInd),'paired');
        % text(maxValueX-(maxValueX-minValueX)*(1/3), (cutoff + 0.75),sprintf('d\'': %.2g',dperiph),...
        %     'Color','w','FontWeight','Bold','FontSize',fontsize)
        % text(maxValueX-(maxValueX-minValueX)*(1/3), (cutoff - 0.75) ,sprintf('d\'': %.2g',dfoveal),...
        %     'Color','w','FontWeight','Bold','FontSize',fontsize)
        % 
        % text(maxValueX-(maxValueX-minValueX)*(1/3), (cutoff - 1.5) ,sprintf('std: %.2g',S),...
        %     'Color','w','FontWeight','Bold','FontSize',fontsize)
        title(sprintf('d\'': %.2g/%.2g(p/f),std:%.2g',dperiph,dfoveal,S))

    case {'x0','y0'}
        line([minValueX, maxValueX], [cutoff,cutoff], 'LineStyle','-', 'Color', 'r', 'LineWidth',1)
        % Here we can do the ecc tests in BarData1 and BarData2
        % Cohens d
        fovealInd = (y <= cutoff);
        dfoveal = computeCohen_d(y(fovealInd), x(fovealInd),'paired');
        dperiph = computeCohen_d(y(~fovealInd), x(~fovealInd),'paired');
        text(maxValueX-(maxValueX-minValueX)*(1/3), (cutoff + 0.75),sprintf('d\'': %.2g',dperiph),...
            'Color','w','FontWeight','Bold','FontSize',fontsize)
        text(maxValueX-(maxValueX-minValueX)*(1/3), (cutoff - 0.75) ,sprintf('d\'': %.2g',dfoveal),...
            'Color','w','FontWeight','Bold','FontSize',fontsize)
    case {'co'}
        % Here we can do the ecc tests in BarData1 and BarData2
        % Cohens d
        d = computeCohen_d(y, x,'paired');
        text(maxValueX-(maxValueX-minValueX)*(1/3), minValueY+fov ,sprintf('d\'': %.2g',d),...
            'Color','w','FontWeight','Bold','FontSize',fontsize)
    otherwise
        error('fieldName not known, only ecc and co implemented ')
        
end





% color things         
% matlab has funky behavior where the size of this influences the size of 
% all future colorbars...
% colormap(zeros(64,3)); 
% cmapValuesHist = colormap('pink');


% caxis([0 maxZ]);
% caxis([0 maxZ/5])
        
end