function  crCreateLinePlot(  R, C_data, cr, ...
                                list_subInds,...
                                list_roiNames,...
                                list_rmNames,...
                                list_rmDescripts,...
                                fname)
% Line plots: from checquerboard to word
% Combine code in figScript_centers_* to create the plots Brian nd I
% discussed
%   - left VOTC, line plot going from word to cb
%   - black if both word and cb are in the same quadrant (); else gray
%   - continuous line if word is inside fovea (0-x (x=2deg)), and CB is out; else dashed
%   - count percentage of continuous black lines over rest:
%      - bin by size (small sizes very unreliable)
%      - 
numRois = length(list_roiNames);
% PLOT IT
for fov      = [ 1.5]
xx = mrvNewGraphWin('LineRadiality and Scatterplot','wide');
set(xx,'Position',[0.005 0.062 .95 .55 ]);
for jj = 1:numRois  
    % data
    ldata = R.L_data{jj}; 
    
    X1 = R.X_rm1{jj};
    Y1 = R.Y_rm1{jj};
    
    X2 = R.X_rm2{jj};
    Y2 = R.Y_rm2{jj};
    
    C = C_data{jj};
    
    ecc_rm1 = R.Ecc_rm1{jj}; 
    ecc_rm2 = R.Ecc_rm2{jj};
    
    ph_rm1 = R.Ph_rm1{jj}; 
    ph_rm2 = R.Ph_rm2{jj};
    
    sm_rm1 = R.Sm_rm1{jj}; 
    sm_rm2 = R.Sm_rm2{jj};
    
    roiName = list_roiNames{jj};
    
  
    % plot on polar map
    %   + black if both word and cb are in the same quadrant (); else gray
    %   + continuous line if word is inside fovea (0-x (x=2deg)), and CB is out; else dashed
    %   + count percentage of continuous black lines over rest:
    %      + bin by size (small sizes very unreliable)
        
    % Filter by size, always, +-inf for no filter
    % Moved this threshold above, now I can't count, make reports in the f
    % sizeMIN  = .5;
    % sizeMAX  = 8;
    
    % sIndw = (sm_rm1 <= sizeMAX) & (sm_rm1 >= sizeMIN); 
    % sIndc = (sm_rm2 <= sizeMAX) & (sm_rm2 >= sizeMIN); 
    % sInd  = (sIndw & sIndc);
    
    % ldata = ldata(sInd); 
    % X1 = X1(sInd);
%     Y1 = Y1(sInd);
%     X2 = X2(sInd);
%     Y2 = Y2(sInd);
%     C = C(sInd);
%     ecc_rm1 = ecc_rm1(sInd); 
%     ecc_rm2 = ecc_rm2(sInd);
%     ph_rm1  = ph_rm1(sInd); 
%     ph_rm2  = ph_rm2(sInd);
%     sm_rm1  = sm_rm1(sInd); 
%     sm_rm2  = sm_rm2(sInd);
%     % histogram(sm_rm1);hold on;histogram(sm_rm2);legend()
    % Obtain angle and ecc again, just in case
    % [PW,ECCW] = cart2pol(X1,Y1);
    % [PC,ECCC] = cart2pol(X2,Y2);
    
   if 0
        fovs   = [0,0.5,1,1.5, 2,3,4,5];
        for ff = 1:length(fovs)
            fov = fovs(ff);
            counter = 0;
            for pp = 1:length(X1)
                pw   = rad2deg(PW(pp));
                pc   = rad2deg(PC(pp));
                qw   = floor(pw/90)+1;
                qc   = floor(pc/90)+1;
                eccw = ECCW(pp);
                eccc = ECCC(pp);

                % if qw==qc && eccw < fov && eccc > fov; counter=counter+1; end
                % if qw==qc && (eccc - eccw) > fov; counter=counter+1; end
                % if qw==qc && (eccc > eccw); counter=counter+1; end
                if (-eccc + eccw) > fov; counter=counter+1; end
                % if eccw < fov && eccc > fov; counter=counter+1; end
                % if qw==qc && ((eccc - eccw) > fov) ; counter=counter+1; end
            end
            if ff==1
                fprintf('\n\n--------------------------------------------------\n');
                % fprintf('SizeMIN: %g, sizeMAX: %g', sizeMIN, sizeMAX);
                % fprintf('  (orig: %g, filtered: %g)', length(sInd), sum(sInd));
                fprintf('\n---------------------------------------------------\n');
            end
            fprintf('Same quadrant and word inside %g deg and cb out, %i out of %i, %0.2f%%\n',...
                    fov,counter,length(X1),counter*100/length(X1))
        end
   end
    testing = false;
    radialityCone=true;
   if 1
        % initialize polar plot
        % Limit plot to visual field circle
        subplot(1,numRois,jj); 
        axis([-cr.defaults.covfig.vfc.fieldRange cr.defaults.covfig.vfc.fieldRange ...
              -cr.defaults.covfig.vfc.fieldRange cr.defaults.covfig.vfc.fieldRange])

        % polar plot
        ff_polarPlot(cr.defaults.covfig.vfc); 
        hold on; 

        % colorbar
        % c = colorbar;
        % colormap(cr.defaults.covfig.vfc.cmapValues)
        % set(c, 'Color', [1 1 1])
        % caxis(cr.defaults.covfig.vfc.cmapRange)
        counter  = 0;
        
        maxang15 = [];
        minang15 = [];
        maxang5  = [];
        minang5  = [];
        for pp = 1:length(X1)
            % pp = 555;
            % lineColor = C(pp,:);
            % pw   = rad2deg(PW(pp));
            % pc   = rad2deg(PC(pp));
            % qw   = floor(pw/90)+1;
            % qc   = floor(pc/90)+1;
            % eccw = ECCW(pp);
            % eccc = ECCC(pp);
            
            pw   = ph_rm1(pp);
            pc   = ph_rm2(pp);
            eccw = ecc_rm1(pp);
            eccc = ecc_rm2(pp);
            
            % [en,an,dn] = crMove2Horiz(e,a,d)

            % if qw==qc ; lineColor = [1 1 1];
            % else; lineColor = [.5 .5 .5]; end

            % if eccw < 5 && eccc > 5 ; lineStyle = '-';
            % else; lineStyle = ':'; end

            % if qw==qc && eccw < 5 && eccc > 5; counter=counter+1; end
            if (abs(eccc - eccw)  > fov)
                if  ((eccc - eccw) > 0) % qw==qc &&
                    counter=counter+1; 
                    lineColor = [.1 .1 .1];
                    lineStyle = ['-'];
                else
                    lineColor = [.5 .5 .5];
                    lineStyle = ['-'];
                end
                
                if radialityCone
                    % WORDS
                    x1 = X1(pp);
                    y1 = Y1(pp);
                    % CB 
                    x2 = X2(pp);
                    y2 = Y2(pp);
                    % We know the polar angle and the ecc, rotate & translate
                    % Rotate angle
                    pcn    = 0;
                    pwn    = pw - pc;
                    % Convert to cartesian
                    [x1n, y1n] = pol2cart(pwn, eccw);
                    [x2n, y2n] = pol2cart(pcn, eccc);
                    % Do the translation only in the x axis
                    if eccc > 5
                        x2nt  = 15;
                        x2dif = x2nt - x2n;
                        x1nt  = x1n + x2dif;
                        % Plot with the same y
                        line([x1nt x2nt], [y1n, y2n], ...
                            'Color', lineColor, 'LineStyle', lineStyle, 'LineWidth',1); 
                    else
                         x2nt  = 5;
                         x2dif = x2nt - x2n;
                         x1nt  = x1n + x2dif;
%                         % Plot with the same y
%                         line([x1nt x2nt], [y1n, y2n], ...
%                             'Color', lineColor, 'LineStyle', lineStyle, 'LineWidth',1); 
                    end
                    
                     % Calculate angles
                     angle = atand(y1n/(x2nt-x1nt));
                    if  ((eccc - eccw) > fov)
                        if eccc > 5
                            if angle > 0
                                maxang15 = [maxang15 angle];
                            end
                            if angle < 0
                                minang15 = [minang15 angle];
                            end
                        else
                            if angle > 0
                                maxang5 = [maxang5 angle];
                            end
                            if angle < 0
                                minang5 = [minang5 angle];
                            end
                        end
                    end
                    
                    if testing
                        line([x1 x2], [y1, y2], ...
                            'Color', 'r', 'LineStyle', lineStyle, ...
                             'LineWidth',1.5); 
                         % Assertions that the calculations are right
                        d      = sqrt((x2-x1)^2 +  (y2-y1)^2)
                        dn     = sqrt((x2nt-x1nt)^2 +  (y2n-y1n)^2)
                    end
                    
                else
                    plot([X1(pp) X2(pp)], [Y1(pp), Y2(pp)], ...
                        'Color', lineColor, 'LineStyle', lineStyle, ...
                         'LineWidth',1.5);       
                end
            end
        end

        % titleName = {
            % ['pRF Center radiality in ' roiName]
            % ['white(' sprintf('%0.2f%% ',counter*100/length(X1)) ...
            %  ')=same quadrant, ecc: (' rmDescript2 ' - ' rmDescript1 ') > ' num2str(fov)]
            % ['sMIN:' num2str(sizeMIN) ',sMAX:' num2str(sizeMAX) ...
            % ' (pRFs ' num2str(sum(sInd)) ' of '  num2str(length(sInd)) ')']
        %     };
        titleName = {
            [strrep(ff_stringRemove(roiName, 'WangAtlas_'),'_','\_')]; 
            % ['slope: ' num2str(meanSlope)];
            % ['ci: ' num2str(ci')];
            % [num2str(numVoxels) ' voxels']
            };
        ttext = strrep(titleName{1},'\_left','');
        
        % title(titleName, 'fontweight', 'bold', 'color', [.1 .1 .1], 'fontsize', 14);
        text(-9, 11, ttext,'Color','k','FontWeight','Bold','FontSize',16)
        
        if ~isempty(maxang15)
            text(-9, 3,sprintf('+%1.1f°',max(maxang15)),'Color','k','FontSize',10)
            text(-9, 1,sprintf('+%1.1f°',median(maxang15)),'Color','r','FontSize',10)
            line([15+5*cosd(180-median(maxang15)),15],[5*sind(median(maxang15)),0],'Color','r')
        end
        if ~isempty(minang15)
            text(-9, -1,sprintf('%1.1f°',median(minang15)),'Color','r','FontSize',10)
            text(-9, -3,sprintf('%1.1f°',min(minang15)),'Color','k','FontSize',10)
            line([15+5*cosd(180-median(minang15)),15],[5*sind(median(minang15)),0],'Color','r')
        end
        
%         if ~isempty(maxang5)
%             text(-9, 3,sprintf('+%1.1f°',max(maxang5)),'Color','k','FontSize',10)
%             text(-9, 1,sprintf('+%1.1f°',median(maxang5)),'Color','r','FontSize',10)
%             line([5+5*cosd(180-median(maxang5)),5],[5*sind(median(maxang5)),0],'Color','r')
%         end
%         if ~isempty(minang5)
%             text(-9, -1,sprintf('%1.1f°',median(minang5)),'Color','r','FontSize',10)
%             text(-9, -3,sprintf('%1.1f°',min(minang5)),'Color','k','FontSize',10)
%             line([5+5*cosd(180-median(minang5)),5],[5*sind(median(minang5)),0],'Color','r')
%         end
        xlim([-15, 18])
        % titlefile = strrep(titleName{1},' ','_');
        % saveas(gcf, fullfile(crRootPath,'local','png',[titlefile '.png']), 'png') 
%         subplot(2,numRois,numRois+jj); 
%         polarhistogram(deg2rad([maxang5 minang5]),100, 'Normalization','probability',...
%             'EdgeAlpha',1, 'EdgeColor','k','FaceAlpha',1, 'FaceColor','k'); 
%         hold on;
%         polarhistogram(deg2rad([maxang15 minang15]),100, 'Normalization','probability',...
%             'EdgeAlpha',1, 'EdgeColor','b','FaceAlpha',.4, 'FaceColor','b'); 
%         legend({'< 5°','> 5°'})
        

   end
   
end

    set(gcf,'color','w');
    if ~isempty(fname)
        saveas(gcf, fullfile(cr.dirs.FIGPNG, [fname '.png']), 'png')
        saveas(gcf, fullfile(cr.dirs.FIGSVG,[fname '.svg']), 'svg')
    end


end

end

