%% plot behavioral scores (towre) against some subject characteristic (coverage extent)
% some hard coding at the moment, will make more genearalizable soon

clear all; close all; clc; 
bookKeeping;

%% modify here

% subjects to include in this plot
list_subInds = [3,5,7:11];

% name of the roi
roiName = 'lVOTRC'; % left_VWFA_rl

% want some ret model information -- dt and model name
dtName = 'Words';
rmName = 'retModel-Words-css.mat';

% name of the behavioral data file. this value will not change
% pathBehavioralData = '/sni-storage/wandell/data/reading_prf/forAnalysis/behavioral/TOWRE_percentileRank.txt';
% pathBehavioralData = '/sni-storage/wandell/data/reading_prf/forAnalysis/behavioral/TOWRE_SWE.txt';
pathBehavioralData = '/sni-storage/wandell/data/reading_prf/forAnalysis/behavioral/TOWRE_PDE.txt';



% coverage contour level, so as to calculate the extent
contourLevel = 0.5; 

% parameters to plot vfc
vfc = ff_vfcDefault; 

%% end modification section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define and intialize some things

numSubs = length(list_subInds);
RFCovAll = zeros(vfc.nSamples, vfc.nSamples, numSubs); 
areaAll = zeros(1, numSubs); 

%% Read in the behavioral data -------------------------------------------
% testing right now. Soon will have to account for the fact that multiple
% tests are stored in this text file
fid = fopen(pathBehavioralData);
S = textscan(fid, '%s%s');  % TODO: get rid of this hard-coding
fclose(fid);

% vector corresponding to subjects of interest
tem = S{2}(list_subInds+1);
dataBehavioral = str2double(tem)
nameBehavioral = S{2}{1};


%% Get the brain data ----------------------------------------------------

%% get the rmroi struct for the subjects we're interested in
RMROI = ff_rmroiMake(roiName, list_subInds, dtName, rmName);

%% get the RFCov matrix from each subject and store it
for ii = 1:numSubs
    
    % this rmroi
    rmroi = RMROI{ii};
    
    % plot the coverage from the individual's rmroi struct
    % [RFcov, figHandle, all_models, weight, data] = rmPlotCoveragefromROImatfile(rm,vfc)
    [RFCov, ~, ~, ~, ~] = rmPlotCoveragefromROImatfile(rmroi,vfc); 
    close; 
    
    RFCovAll(:,:,ii) = RFCov; 
    
end

%% get the individual coverage area and store it

for ii = 1:numSubs
    % this subject's coverage
    RFCov = RFCovAll(:,:,ii); 
    
    [covAreaPercent, covArea] = ff_coverageArea(contourLevel, vfc, RFCov); 
    areaAll(ii) = covArea;
    
end



%% plot 
close all; 
figure; 
hold on; 

for ii = 1:numSubs
    
    subInd = list_subInds(ii);
    subColor = list_colorsPerSub(subInd, :);
    plot(areaAll(ii), dataBehavioral(ii),'o', 'MarkerFaceColor', subColor, ...
        'MarkerEdgeColor', [0 0 0], 'MarkerSize', 18, 'LineWidth',2);
    
end

%% Calculate the correlation and parameters of the line
rho = corrcoef(areaAll, dataBehavioral);
corrValue = rho(1,2);

p = polyfit(areaAll', dataBehavioral,1);
xlim = get(gca, 'Xlim'); 
lin = linspace(xlim(1),xlim(2)); 
plot(lin, polyval(p,lin), 'LineWidth', 5, 'Color', [.1 .1 .1])

%% Plot properties
grid on; 
xlabel(['Cov area (deg^2) (' num2str(contourLevel) ' contour level)'])
ylabel(nameBehavioral)
titleName = {
    [nameBehavioral ' and coverage area']
    ['rho: ' num2str(corrValue)]
    };
title(titleName, 'fontweight', 'bold')
ff_dropboxSave;

ff_legendSubject(list_subInds);
title('legend')
ff_dropboxSave; 

