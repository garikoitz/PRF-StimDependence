function  [percentAboveSubs,perROI] = crCreateScatterplot(R, C_data, cr, ...
                                                          list_subInds,...
                                                          list_roiNames,...
                                                          list_rmDescripts,...
                                                          fieldName, ...
                                                          fname)
%%
% colormap for histogram
% cmapValuesHist = colormap('pink');
% cmapValuesHist_tem = colormap('hot');
% cmapValuesHist = cmapValuesHist_tem(2:55, :); 
colormap(zeros(64,3)); % matlab has funky behavior where the size of this influences the size of all future colorbars...
cmapValuesHist = colormap('pink');
close; 

% whether looking at a subject by subject basis
% subIndividually = false; 

% field to plot. Ex:  
% 'co': variance explained 
% 'ecc': eccentricity 
% 'sigma': effective size 
% 'sigma1': sigma major
% 'numvoxels' for number of voxels in roi
% fieldToPlotDescript is for axis labels and plot titles
%     'sigma'       : effective sigma
%     'sigma1'      : sigma major
%     'ecc'         : eccentricity
%     'co'          : variance explained 
%     'exponent'    : exponent
%     'betaScale'   : how much to scale the predicted tseries by
%     'meanMax'     : mean of the top 8 values
%     'meanPeaks'   : mean of the outputs of matlab's meanPeaks
list_fieldNames  = {    
    'co'
    'ecc'
%     'sigma'
%     'ph'
    }; 

list_fieldDescripts = {
    'variance explained'
    'eccentricity'
%     'sigma';afo
%     'polar angle'
    }; 

% which plots do we want? lots we can make ...
plot_fit = false; % plotting the across-subject bootstrapped line w/ CIs

% transparency of the plots
alphaValue = 0.4; 



% location of the colorbar
% default: 'eastoutside'
% 'southoutside': 
cbarLocation = 'eastoutside';

% end modification section

numSubs = length(list_subInds);
numRois = length(list_roiNames);

% number of fields
% numFields = length(list_fieldNames);

% rm descriptions
rm1Descript = list_rmDescripts{1}; 
rm2Descript = list_rmDescripts{2}; 

% initialize structs / matrices for mixed effects
% subjectLines = cell(numSubs, numRois, numFields); % because pairwise
subjectLines = cell(numSubs, numRois); % because pairwise

% initialize struct for calculating the percentage of voxels above the
% identity line
% percentAbovePooled = zeros(numRois, numFields);
percentAbovePooled = zeros(numRois);
% percentAboveSubs   = zeros(numSubs, numRois, numFields);
percentAboveSubs   = zeros(numSubs, numRois);
diffAboveSubs      = zeros(numSubs, numRois);
% A = cell(numFields*numRois, 5);
A = cell(numRois, 5);

% get the cell of rms so that we can threshold
% rmroiCell = ff_rmroiCell(cr, list_subInds, list_roiNames, list_dtNames, list_rmNames, ...
%     'list_path', list_path);

% Threshold and get identical voxels for each subject
% In comparing ret models, the collection of voxels may not be the same
% because of the thresholding. In this cell we redefine the rmroi
% R.rmroiCellSameVox = cell(size(rmroiCell));

% for jj = 1:numRois
%     for ii = 1:numSubs        
        % get identical voxels for each subject's roi over all ret models
%         D = rmroiCell(ii,jj,:);
%         R.rmroiCellSameVox(ii,jj,:) = ff_rmroiGetSameVoxels(D, vfc);        
%     end
% end

% close all;


    % subplot(2,numRois,numRois+jj); 

    
    % field-specific properties
    % fieldName = list_fieldNames{ff}; % We pass it now...
    % fieldNameDescript = list_fieldDescripts{ff}; 
    radius = 1;
    if strcmp(fieldName, 'sigma1') 
        maxValue = cr.defaults.covfig.vfc.sigmaMajthresh(2);
    elseif strcmp(fieldName, 'sigma')          
        maxValue = cr.defaults.covfig.vfc.sigmaEffthresh(2);
    elseif strcmp(fieldName, 'ecc')
        maxValue = cr.defaults.covfig.vfc.eccthresh(2);
        fov = 1.5; % width of the band
        nrows = 2; ncols = 3;
        position = [0.005 0.062 .95 .7 ];
        radius = 2;
    elseif strcmp(fieldName, 'co')
        maxValue = 1; 
        fov = 0.2; % width of the band
        % nrows = 1; ncols = 6;
        nrows = 2; ncols = 3;
        % position = [0.005 0.062 .95 .6 ];
        position = [0.005 0.062 .95 .7 ];
        radius = .15;
    elseif strcmp(fieldName, 'exponent')
        maxValue = 2; 
    elseif strcmp(fieldName, 'meanMax')
        maxValue = 20; 
    elseif strcmp(fieldName, 'meanPeaks')
        maxValue = 10; 
    elseif strcmp(fieldName, 'betaScale')
        maxValue = 5; 
    elseif strcmp(fieldName, 'x0') || strcmp(fieldName, 'y0')
        maxValue = cr.defaults.covfig.vfc.fieldRange;
    elseif strcmp(fieldName, 'ph')
        maxValue = 2*pi; 
    else
        error('Define the maxValue so we can normalize and fit the beta distribution.');
    end


    xx = mrvNewGraphWin('Scatterplots');
    set(xx,'Position',position);
    for jj = 1:numRois
        roiName = list_roiNames{jj};
        subplot(nrows,ncols,jj);
    

        
        axisLims = [0 maxValue]; 
        BarData1 = [];
        BarData2 = [];
        
        poolallsubs1=[];
        poolallsubs2=[];
        for ii = 1:numSubs

            subInd = list_subInds(ii);
            
            % rmRois for different ret models
            rmroi1 = R.rmroiCellSameVox{ii,jj,1}; 
            rmroi2 = R.rmroiCellSameVox{ii,jj,2};

             if ~isempty(rmroi1)
                 
                % get the data
                x1 = eval(['rmroi1.' fieldName]);
                x2 = eval(['rmroi2.' fieldName]);

                % shift so that the smallest value is 0
                if strcmp(fieldName, 'ph')
                    x1 = x1 + pi; 
                    x2 = x2 + pi;
                end
                
                %% For the mixed effects
                % fit a line for each subject
                % p = polyfit(x1, x2, 1);              
                % subjectLines{ii,jj,ff} = p;                 
                % [B, BINT] = regress(x2', x1'); 
                % b.slope     = B; 
                % b.ci95      = BINT; 
                
                % p = polyfit(x1, x2, 1);
                % b.pintercept = p(1);
                % b.pslope = p(2); 
                
                % subjectLines{ii,jj,ff} = b; 
                
                
                poolallsubs1=[poolallsubs1, x1];
                poolallsubs2=[poolallsubs2, x2];
                
                
                %% the percentage of voxels above the identityLine
                perAbove = sum(x2 > x1) / length(x2);
                % percentAboveSubs(ii,jj,ff) = perAbove; 
                percentAboveSubs(ii,jj) = perAbove; 
                
                % Just store the data and put it outside
%                 diffAboveSubs(ii,jj) = x2-x1; 
                
                %% concatenate
                BarData1 = [BarData1, x1];
                BarData2 = [BarData2, x2];
           
            end
        end % end loop over subjects
        
        %% mixed effects: fit a line to individual subjects
        % slopes = nan(1, numSubs);   
        % slopesp = nan(1,numSubs);
        % interceptsp = nan(1,numSubs);
        % percents = percentAboveSubs(:,jj,ff);
        
        % for ii = 1:numSubs
        %     b = subjectLines{ii,jj,ff}; 
        %     if ~isempty(b)
        %         slopes(ii) = b.slope; 
        %     end
        % end

        % the calculating. nan will cause bootci to error
        % slopes(isnan(slopes)) = []; 
        % percents(isnan(percents)) = [];
        
        % table things. 
        % (1)roiName (2)fieldName (3) ciLow (4)ciHigh (5)mean
        % tind = (jj-1)*numFields + ff; 
        
        
        % if numSubs > 1
            % numbs = 1000; 
            % [ci, bootstat] = bootci(numbs, @mean, slopes);
            % meanSlope = mean(bootstat); 
            
            % % table things for percent above identityLine
            % [ciPer, bootstatPer] = bootci(numbs, @mean, percents);
            % A{tind,1} = ff_stringRemove(roiName, '_rl');
            % A{tind,2} = fieldName;
            % A{tind,3} = ciPer(1);
            % A{tind,4} = ciPer(2);
            % A{tind,5} = mean(bootstatPer);
            
        % else
        %     meanSlope = nan; 
        %     ci = nan; 
        %   end
                
        %% calculations related to the percentage above the identityLine
        % percentabove = sum(BarData2 > BarData1) / length(BarData1);
        % percentAbovePooled(jj,ff) = percentabove;
        
        % properties related to both types of scatter plots
        % coloring by number of voxels or percentage of voxels
        npoints = 100; 
        
        % 3d histogram heat map -- absolute number of voxels
        ff_histogramHeat(BarData1, BarData2, maxValue, maxValue,radius,cmapValuesHist,fov,roiName);
        numVoxels = length(BarData1); 
        % axes and title
        switch fieldName
            case {'ecc'}
                if jj==1 || jj==4
                    ylabel(['pRF eccentricity for ' rm2Descript ' (deg)'])
                end
                if jj>3
                    xlabel(['pRF eccentricity for ' rm1Descript ' (deg)'])
                end
            case {'co'}
                if jj==1 || jj==4
                    ylabel(['Variance explained for ' rm2Descript ' (%)'])
                end
                if jj>3
                    xlabel(['Variance explained for ' rm1Descript ' (%)'])
                end
            otherwise
                if jj==1 || jj==4
                    ylabel(['' rm2Descript ''])
                end
                if jj>3
                    xlabel(['' rm1Descript ''])
                end
        end
        
    end % loop over rois

    % titleName = {
    %     [strrep(ff_stringRemove(roiName, 'WangAtlas_'),'_','\_') '.' fieldName];
    % ['slope: ' num2str(meanSlope)];
    % ['ci: ' num2str(ci')];
    % [num2str(numVoxels) ' voxels']
    %      };
    % title(titleName, 'FontWeight', 'Bold');
    % fname = [titleName{1} '_' fieldName '_band-2x' num2str(fov)];
    if ~isempty(fname)
        set(gcf, 'InvertHardcopy', 'off')
        saveas(gcf, fullfile(cr.dirs.FIGPNG, [fname '.png']), 'png')
        saveas(gcf, fullfile(cr.dirs.FIGSVG,[fname '.svg']), 'svg')
    end

end

