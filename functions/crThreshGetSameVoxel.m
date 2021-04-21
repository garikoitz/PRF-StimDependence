function [rmroiCellSameVox,C_data,cmapValuesHist,maxValue]=crThreshGetSameVoxel(...
                                                                cr,...
                                                                rmroiCell, ...
                                                                list_subInds,...
                                                                list_roiNames,...
                                                                list_rmNames,...
                                                                list_rmDescripts,...
                                                                fname)

%  Threshold and get identical voxels for each subject
% values to threshold the RM struct by
% vfc = ff_vfcDefault_Hebrew;
% The defaults are different for the two projects, so for now use the most
% restrictive one
% vfc threshold
vfc.prf_size         = true;
vfc.fieldRange       = 15; % 7;
vfc.method           = 'max';
vfc.newfig           = true;
vfc.nboot            = 50;
vfc.normalizeRange   = true;
vfc.smoothSigma      = false;
% Thresholds
vfc.sigmaEffthresh   = [.2 15];  % [0  7]; % sigma effect (sigmaMajor/sqrt(exponent))
vfc.sigmaMajthresh   = [0.5 8];   % [0 14]; % sigma major (before the exponent)
vfc.cothresh         = 0.20; % Variance explained, 20%
vfc.cothreshceil     = 1; % 0.2 looking for noise. 1 for normal thresholding. don't get voxels higher than this
vfc.threshByCoh      = true;
vfc.eccthresh        = [0 vfc.fieldRange];
vfc.quadthresh       = [1 4];  % If center not on these quads, remove
% Fig control
vfc.nSamples         = 128;
vfc.meanThresh       = 0;
vfc.weight           = 'fixed';
vfc.weightBeta       = false;
vfc.cmap             = 'jet';
vfc.clipn            = 'fixed';
vfc.addCenters       = false;
vfc.verbose          = prefsVerboseCheck;
vfc.dualVEthresh     = 0;
vfc.ellipsePlot      = false;
vfc.ellipseLevel     = 0.5;
vfc.ellipseColor     = [1 0 0];
vfc.contourPlot      = true;
vfc.contourLevel     = 0.5;
vfc.contourColor     = [0 0 0];
vfc.tickLabel        = false;
vfc.contourBootstrap = false;
vfc.gridColor        = [.6 .6 .6];
vfc.backgroundColor  = [1 1 1];
vfc.cmapRange        = [0 pi]; % the range over which color bar
vfc.cmapValues       = flipud(jetCmap(0,128)); % the colorbar values
vfc.alphaValue       = ''; % 0.5;
vfc.alphaValueDot    = ''; % 0.8;
vfc.lineWidth        = 1.5; % thicker lines for transparent. 1 works well for opaque
% by definition, when eccentricity does not shift a lot, theta will not
% shift either. We can look at theta shifts in the voxels whose
% eccentricity have shifted a lot
% thetaShiftByEccThresh = true;
% eccThresh = 3;
% Update the cr structure that I created that go into the new functions
cr.defaults.covfig.vfc = vfc;
% INITIALIZE SOME THINGS
numRois = length(list_roiNames);
numSubs = length(list_subInds);
% cell for linearizing the data (a vector for each ROI)
L_data  = cell(1, numRois);
X_rm1   = cell(1, numRois);
Y_rm1   = cell(1, numRois);
X_rm2   = cell(1, numRois);
Y_rm2   = cell(1, numRois);
Ecc_rm1 = cell(1, numRois);
Ecc_rm2 = cell(1, numRois);
rmDescript1 = list_rmDescripts{1};
rmDescript2 = list_rmDescripts{2};
% In comparing ret models, the collection of voxels may not be the same
% because of the thresholding. In this cell we redefine the rmroi
rmroiCellSameVox = cell(size(rmroiCell));
for jj = 1:numRois
    for ii = 1:numSubs        
        % get identical voxels for each subject's roi over all ret models
        D = rmroiCell(ii,jj,:);
        % GLU EDIT function: remove voxels from the oppossite hemifield
        rmroiCellSameVox(ii,jj,:) = ff_rmroiGetSameVoxels(D, cr.defaults.covfig.vfc);
    end
end

% Linearize the data
% Take the difference between 2 rms. 
% Also store the x and y data
for jj = 1:numRois
    % initializing the difference of the centers' thetas
    ldata = []; 

    % intializing the location of the centers
    xdata_rm1 = [];
    ydata_rm1 = []; 
    xdata_rm2 = [];
    ydata_rm2 = [];
    
    % initializing eccentrcity
    ecc_rm1   = [];
    ecc_rm2   = []; 
    
    % initializing angle
    ph_rm1    = [];
    ph_rm2    = []; 
    
    % initializing size
    sm_rm1    = [];
    sm_rm2    = []; 
    
    for ii = 1:numSubs
        rmroi1 = rmroiCellSameVox{ii,jj,1};
        rmroi2 = rmroiCellSameVox{ii,jj,2};
        
        % some subjects don't have 
        if ~isempty(rmroi1) & ~isempty(rmroi2)
            data1 = rmroi1.ph;
            data2 = rmroi2.ph;

            % the difference between centers' thetas.
            % this will determine the color of the line
            % we take absolute value because we are interested in the magnitude
            % of the rotation and not the direction
            fieldDiffOver = abs(data2 - data1);  

            % Note that the difference will range between 0 and 2pi. 
            % We want to constrain values to be between and pi (again not 
            % interested in the direction of the rotation but the magnitude of it)
            % For values greater than pi, subtract it from 2pi
            fieldDiff = ff_polarAngleBetween0AndPi(fieldDiffOver);

            ldata = [ldata fieldDiff];

            % the location of the pRF centers
            xdata_rm1 = [xdata_rm1 rmroi1.x0]; 
            ydata_rm1 = [ydata_rm1 rmroi1.y0]; 

            xdata_rm2 = [xdata_rm2 rmroi2.x0]; 
            ydata_rm2 = [ydata_rm2 rmroi2.y0]; 

            ecc_rm1   = [ecc_rm1 rmroi1.ecc];
            ecc_rm2   = [ecc_rm2 rmroi2.ecc];
            
            ph_rm1    = [ph_rm1 rmroi1.ph];
            ph_rm2    = [ph_rm2 rmroi2.ph];
            
            sm_rm1    = [sm_rm1 rmroi1.sigma1];
            sm_rm2    = [sm_rm2 rmroi2.sigma1];
        end        
    end
    L_data{jj} = ldata; 
    
    X_rm1{jj}   = xdata_rm1;
    Y_rm1{jj}   = ydata_rm1;
    
    X_rm2{jj}   = xdata_rm2;
    Y_rm2{jj}   = ydata_rm2;
    
    Ecc_rm1{jj} = ecc_rm1; 
    Ecc_rm2{jj} = ecc_rm2;
    
    Ph_rm1{jj}  = ph_rm1; 
    Ph_rm2{jj}  = ph_rm2;
    
    Sm_rm1{jj}  = sm_rm1; 
    Sm_rm2{jj}  = sm_rm2;
    
end
% Get a colormap according to the linearized data in L_data
for jj = 1:numRois    
    ldata = L_data{jj};     
    cdata = ff_colormapForValues(ldata, cr.defaults.covfig.vfc.cmapValues, ...
                                        cr.defaults.covfig.vfc.cmapRange);    
    C_data{jj} = cdata;  
end

% colormap for histogram
% cmapValuesHist = colormap('pink');
% cmapValuesHist_tem = colormap('hot');
% cmapValuesHist = cmapValuesHist_tem(2:55, :); 
colormap(zeros(64,3)); % matlab has funky behavior where the size of this influences the size of all future colorbars...
cmapValuesHist = colormap('pink');
close; 

% whether looking at a subject by subject basis
subIndividually = false; 

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
numRms  = length(list_rmNames);

% number of fields
numFields = length(list_fieldNames);

% rm descriptions
rm1Descript = list_rmDescripts{1}; 
rm2Descript = list_rmDescripts{2}; 

% initialize structs / matrices for mixed effects
subjectLines = cell(numSubs, numRois, numFields); % because pairwise

% initialize struct for calculating the percentage of voxels above the
% identity line
percentAbovePooled = zeros(numRois, numFields);
percentAboveSubs   = zeros(numSubs, numRois, numFields);
A = cell(numFields*numRois, 5);

% get the cell of rms so that we can threshold
% rmroiCell = ff_rmroiCell(cr, list_subInds, list_roiNames, list_dtNames, list_rmNames, ...
%     'list_path', list_path);

% Threshold and get identical voxels for each subject
% In comparing ret models, the collection of voxels may not be the same
% because of the thresholding. In this cell we redefine the rmroi
% rmroiCellSameVox = cell(size(rmroiCell));

% for jj = 1:numRois
%     for ii = 1:numSubs        
        % get identical voxels for each subject's roi over all ret models
%         D = rmroiCell(ii,jj,:);
%         rmroiCellSameVox(ii,jj,:) = ff_rmroiGetSameVoxels(D, vfc);        
%     end
% end

% close all;


    % subplot(2,numRois,numRois+jj); 

    
for ff = 1:numFields

    % field-specific properties
    fieldName = list_fieldNames{ff};
    fieldNameDescript = list_fieldDescripts{ff}; 

    if strcmp(fieldName, 'sigma1') 
        maxValue = cr.defaults.covfig.vfc.sigmaMajthresh(2);
    elseif strcmp(fieldName, 'sigma')          
        maxValue = cr.defaults.covfig.vfc.sigmaEffthresh(2);
    elseif strcmp(fieldName, 'ecc')
        maxValue = cr.defaults.covfig.vfc.eccthresh(2);
        fov = 1.5; % width of the band
        nrows = 2; ncols = 3;
        position = [0.005 0.062 .9 .7 ];
    elseif strcmp(fieldName, 'co')
        maxValue = 1; 
        fov = 0.2; % width of the band
        nrows = 1; ncols = 6;
        position = [0.005 0.062 .9 .5 ];
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

        for ii = 1:numSubs

            subInd = list_subInds(ii);
            
            % rmRois for different ret models
            rmroi1 = rmroiCellSameVox{ii,jj,1}; 
            rmroi2 = rmroiCellSameVox{ii,jj,2};

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
                
                %% the percentage of voxels above the identityLine
                perAbove = sum(x2 > x1) / length(x2);
                percentAboveSubs(ii,jj,ff) = perAbove; 
                
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
        ff_histogramHeat(BarData1, BarData2, maxValue, maxValue, 50,cmapValuesHist,fov,roiName);
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
                if jj==1 % || jj==4
                    ylabel(['Variance explained for ' rm2Descript ' (%)'])
                end
                % if jj>3
                    xlabel(['Variance explained for ' rm1Descript ' (%)'])
                % end
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
    saveas(gcf, fullfile(cr.dirs.FIGPNG, [fname '.png']), 'png')
    saveas(gcf, fullfile(cr.dirs.FIGSVG,[fname '.svg']), 'svg')
end % loop over fields

end
