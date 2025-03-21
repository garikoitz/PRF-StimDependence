%% SOME THINGS TO CHECK HERE
% use summary_behavioral_vsCoveragePlot

%% x axis: behavioral score. y axis: some coverage parameter
close all; clear all; clc; 
bookKeeping; 

%% modify here

markerSize = 18; % 11

% the behavior score to use. options:
% 'towre_percentile', 'towre_swe', 'towre_pde', ...
list_behaviors = {
    'towre_swe'
    'towre_pde'
%     'towre_percentile'
    };

% list of subjects with behavioral score to use
list_subInds = [3,5,7:11];

% the coverage parameter to use. options:
% 'coverageArea'
list_coverageParams = {
    'coverageArea'
%     'coverageAreaPercent'
    };

% -1 if we want to use all values
contourLevel = -1; 

% path of behavior .mat files
list_behaviorData = {
    '/sni-storage/wandell/data/reading_prf/forAnalysis/data/behavioral/towre/towre.mat'
    };

% rois we want to look at
list_roiNames = {
    'lVOTRC'
%     'right_VWFA_rl'
%     'combined_VWFA_rl'
    };

% ret models we want to look at
list_dtNames = {
    'Words';
    };
list_rmNames = {
    'retModel-Words-css.mat'
    };

% vfc threshold -- maybe we should store this somewhere and then load
vfc = ff_vfcDefault; 


%% define things
numSubs = length(list_subInds);
numBehaviors = length(list_behaviors);
numCoverageParams = length(list_coverageParams);

coverageParamMat = zeros(numSubs, numCoverageParams);
behaviorMat = zeros(numSubs, numBehaviors);
list_colors = list_colorsPerSub;

%% do things -----------------

%% load all the behavioral matrices
for bb = 1:length(list_behaviorData)
    load(list_behaviorData{bb})
end

% get the matrix of behavioral scores


for bb = 1:numBehaviors
    behavior = list_behaviors{bb};
    thedata = eval(behavior);

    for ii = 1:numSubs
        subInd = list_subInds(ii);
        behaviorMat(ii,bb) = thedata{subInd,1};
    end    
end

%% get the rmroi cell for the subjects
rmroiCell = ff_rmroiCell(list_subInds, list_roiNames, list_dtNames, list_rmNames);

% and then get a cell of RFcov
rfcovCell = cell(1, numSubs);

for ii = 1:numSubs
   rfcovCell{ii} = rmPlotCoveragefromROImatfile(rmroiCell{ii}, vfc); 
end

%% get the coverageParams

for ii = 1:numSubs
    
    rfcov = rfcovCell{ii};
    
    for cc = 1:numCoverageParams
        coverageParam = list_coverageParams{cc};

        switch coverageParam
            case 'coverageArea'
                [covAreaPercent, covArea] = ff_coverageArea(contourLevel, vfc, rfcov);
                coverageParamMat(ii,cc) = covArea;
            case 'coverageAreaPercent'
                [covAreaPercent, covArea] = ff_coverageArea(contourLevel, vfc, rfcov);
                coverageParamMat(ii,cc) = covAreaPercent;
        end            
    end
end

%% loop over the behavioral tests and make plots

for bb = 1:numBehaviors
    
    behavior = list_behaviors{bb};
    
    for cc = 1:numCoverageParams
        figure; hold on;
        coverageParam = list_coverageParams{cc};
        
        for ii = 1:numSubs
            subInd = list_subInds(ii);
            subColor = list_colors(subInd,:);
            subBehavior = behaviorMat(ii,bb);
            subCoverage = coverageParamMat(ii,cc);
            h(ii) = plot(subCoverage, subBehavior, 'Marker', 'o', ...
                'MarkerSize', markerSize, ...
                'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', subColor, ...
                'LineWidth',2.5, ...
                'Color', 'k');
       
        end
                
        %% plot line
        coverageData = coverageParamMat(:,cc);
        behaviorData = behaviorMat(:,cc);
        
        p = polyfit(coverageData, behaviorData,1); 
        xlim = get(gca, 'Xlim'); 
        lin = linspace(xlim(1),xlim(2)); 
        plot(lin, polyval(p,lin), 'LineWidth', 5, 'Color', [.1 .1 .1])
        
        % the correlation
        rho = corrcoef(coverageData, behaviorData);
        corrValue = rho(1,2);
        
        % plot properties
        L = legend(h, list_sub{list_subInds});
        set(L,'Location', 'NorthEastOutside')
        
        grid on;
        xlabel(coverageParam, 'FontWeight', 'Bold')
        ylabel(behavior, 'FontWeight', 'Bold')
        titleName = {
            [behavior ' vs. ' coverageParam]
            ['rho: ' num2str(corrValue)]
            };
        title(titleName, 'FontWeight', 'Bold', 'FontSize',18)
        ff_dropboxSave;
        
    end
    
end
