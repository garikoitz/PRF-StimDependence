function  [diffs15, diffs5,...
           posangles15,negangles15,...
           posangles5,negangles5] = crCreateLinePlot(R, cr,list_roiNames,fname)

% Line plots: from checquerboard to word
% Combine code in figScript_centers_* to create the plots Brian nd I
% discussed
%   - left VOTC, line plot going from word to cb
%   - black if both word and cb are in the same quadrant (); else gray
%   - continuous line if word is inside fovea (0-x (x=2deg)), and CB is out; else dashed
%   - count percentage of continuous black lines over rest:
%      - bin by size (small sizes very unreliable)
%      - 
numRois     = length(list_roiNames);
diffs15     = cell(numRois,1);
diffs5      = cell(numRois,1);
posangles15 = cell(numRois,1);
negangles15 = cell(numRois,1);
posangles5  = cell(numRois,1);
negangles5  = cell(numRois,1);
% PLOT IT
for fov      = [0]
xx = mrvNewGraphWin('LineRadiality and Scatterplot','wide');
set(xx,'Position',[0.005 0.062 .95 .55 ]);

for necc5check=1:2
    
for jj = 1:numRois  
    % data
    ecc_rm1 = R.Ecc_rm1{jj}; 
    ecc_rm2 = R.Ecc_rm2{jj};
    
    ph_rm1 = R.Ph_rm1{jj}; 
    ph_rm2 = R.Ph_rm2{jj};
    
    % Separate the less than and more than 5 degs, see scattersplots
    if necc5check == 1
        ecc_rm1_sub5ind = (ecc_rm2 <= 5);
        % Filter the rest of the data
        ecc_rm1 = ecc_rm1(ecc_rm1_sub5ind); 
        ecc_rm2 = ecc_rm2(ecc_rm1_sub5ind);
        
        ph_rm1 = ph_rm1(ecc_rm1_sub5ind);
        ph_rm2 = ph_rm2(ecc_rm1_sub5ind);
    else
        ecc_rm1_above5ind = (ecc_rm2 > 5);
        % Filter the rest of the data
        ecc_rm1 = ecc_rm1(ecc_rm1_above5ind); 
        ecc_rm2 = ecc_rm2(ecc_rm1_above5ind);
        
        ph_rm1 = ph_rm1(ecc_rm1_above5ind);
        ph_rm2 = ph_rm2(ecc_rm1_above5ind);
    end
    
    roiName = list_roiNames{jj};
    
    % initialize polar plot
    % Limit plot to visual field circle
    subplot(2,numRois,(necc5check-1)*numRois + jj);
    
    axis([-cr.defaults.covfig.vfc.fieldRange cr.defaults.covfig.vfc.fieldRange ...
        -cr.defaults.covfig.vfc.fieldRange cr.defaults.covfig.vfc.fieldRange])
    
    % polar plot
    % ff_polarPlot(cr.defaults.covfig.vfc);
    radius = 15;
    minAng = -25;
    maxAng = 25;
    lightGray = [.6 .6 .6 ];
    plot([0,radius],[0,0],'Color',lightGray,'LineWidth',2,'LineStyle','-');
    hold on; axis equal
    plot([0,radius*cosd(maxAng)],[0,radius*sind(maxAng)],'Color',lightGray,'LineWidth',2,'LineStyle','-')
    plot([0,radius*cosd(minAng)],[0,radius*sind(minAng)],'Color',lightGray,'LineWidth',2,'LineStyle','-')
    % plot the arc
    angles = minAng:.1:maxAng;
    x = zeros(length(angles),1); y = x;
    for ii = 1:length(angles)
        x(ii) = radius * cosd(angles(ii));
        y(ii) = radius * sind(angles(ii));
    end
    plot (x,y,'Color',lightGray,'LineWidth',2,'LineStyle','-')
    
    % Initialize strcutures per ROI to pass to output
    diffWordCB = [];
    posang     = [];
    negang     = [];
    for pp = 1:length(ph_rm1)
        pw   = ph_rm1(pp);
        pc   = ph_rm2(pp);
        eccw = ecc_rm1(pp);
        eccc = ecc_rm2(pp);

        % fov can be 1.5, to exclude all values inside the green band of the
        % scatterplots
        if (abs(eccc - eccw)  > fov)
            % Rotate angle
            pcn    = 0;
            pwn    = pw - pc;
            % Convert to cartesian
            [x1n, y1n] = pol2cart(pwn, eccw);
            [x2n, y2n] = pol2cart(pcn, eccc);
            % Do the translation only in the x axis
            x2nt  = 15;
            x2dif = x2nt - x2n;
            x1nt  = x1n + x2dif;
                        
            % Store the values
            if (x2nt-x1nt) > 0
            % if (eccc-eccw) > 0
                % Store values
                diffWordCB = [diffWordCB, abs(eccc-eccw)];
                % Line color
                lineColor = [.1 .1 .1];
                lineStyle = ['-'];
                % Calculate angles
                angle = atand(y1n/(x2nt-x1nt));
                
                if angle > 0
                    posang = [posang angle];
                end
                if angle < 0
                    negang = [negang angle];
                end
                
            else
                % Store values
                diffWordCB = [diffWordCB, -abs(eccc-eccw)];
                % Line Color
                lineColor = [.5 .5 .5];
                lineStyle = ['-'];
            end
            % Plot this line
            line([x1nt x2nt], [y1n, y2n], ...
                'Color', lineColor, 'LineStyle', lineStyle, 'LineWidth',.5);
            
        end
    end
    
    ylim([radius*sind(minAng), radius*sind(maxAng)])
    xlim([0, 17.5])
    set(gca,'ycolor','w')
    set(gca,'xcolor','w')

    % Plot ROI in each subplot
    titleName = {[strrep(ff_stringRemove(roiName, 'WangAtlas_'),'_','\_')];};
    ttext     = strrep(titleName{1},'\_left','');
    text(0, 3, ttext,'Color','k','FontWeight','Bold','FontSize',16)
    % Plot median values in each subplot
    if ~isempty(posang)
        % text(5, 3,sprintf('+%1.1f°',max(posang)),'Color','k','FontSize',10)
        text(4, 1,sprintf('+%1.1f°',median(posang)),'Color','r','FontSize',14)
        line([radius+5*cosd(180-median(posang)),15],[5*sind(median(posang)),0],'Color','r')
    end
    if ~isempty(negang)
        text(4, -1,sprintf('%1.1f°',median(negang)),'Color','r','FontSize',14)
        % text(5, -3,sprintf('%1.1f°',min(negang)),'Color','k','FontSize',10)
        line([radius+5*cosd(180-median(negang)),15],[5*sind(median(negang)),0],'Color','r')
    end
    
    if necc5check == 2
        ttext = 'Checkers eccentricity > 5 deg';
        if jj==1
            text(0, 7.5 , ttext,'Color','k','FontWeight','Bold','FontSize',16)
        end
        % Add values per ROI to be used outside this script
        diffs15{jj}=diffWordCB;
        posangles15{jj}=posang;
        negangles15{jj}=negang;
    else
        
        ttext = 'Checkers eccentricity <= 5 deg';
        if jj==1
            text(0, 7.5, ttext,'Color','k','FontWeight','Bold','FontSize',16)
        end
        % Add values per ROI to be used outside this script
        diffs5{jj}=diffWordCB;
        posangles5{jj}=posang;
        negangles5{jj}=negang;
    end
    
end % numRois

end % ecc check
    
    % SAVE THE FIG
    set(gcf,'color','w');
    if ~isempty(fname)
        saveas(gcf, fullfile(cr.dirs.FIGPNG, [fname '.png']), 'png')
        saveas(gcf, fullfile(cr.dirs.FIGSVG,[fname '.svg']), 'svg')
    end

end % fovs-s

end % function

