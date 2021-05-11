function figFunction_coverage_maxProfile_group(cr, subinds, varargin)

% GLU created this function from the script figScript_coverage...

%% Reading circuitry field of view paper (2016). Making figures.
% Uses this script: coverage_plotGroupAverage.m (should make into a
% function)
% (this script used to be called: figScript_coverage_maxProfile_group)

% clear all; close all; clc; 
% bookKeeping;


% Make varargin lower case, remove white spaces...
varargin = mrvParamFormat(varargin);
% Parse
p = inputParser;
p.addRequired('cr'        , @isstruct);
p.addRequired('subinds'   , @isnumeric);

p.addOptional('flip'      , true,  @islogical);
p.addOptional('rmroicell' , {}, @iscell);
p.addOptional('vers','v01', @ischar);
p.addOptional('fname',''  , @ischar);
p.addOptional('invisible',false  , @islogical);


% Parse. Assign result inside each case
p.parse(cr, subinds, varargin{:});
% Read here only the generic ones
flip      = p.Results.flip;
rmroiCell = p.Results.rmroicell;
vers      = p.Results.vers;
fname     = p.Results.fname;
invisible = p.Results.invisible;

if isempty(fname);savefig=false;else;savefig=true;end


%%
if invisible
    visibility = 'off';
else
    visibility = 'on';
end



%%
vfc = cr.defaults.covfig.vfc;
% visual field plotting thresholds
% vfc = ff_vfcDefault; 
% vfc.cothresh = 0.2; 
% vfc.method = 'max'; % avg | 'max'

% session list, see bookKeeping
list_path     = cr.bk.list_sessionRet; 

% subjects to do this for, see bookKeeping
% %[31:36 38 39:44] % Hebrew
list_subInds  = subinds; 

% whether we want to plot the half max contour
plotContour   = vfc.contourPlot; 
contourLevel  = vfc.contourLevel; 

% rois we want to look at
list_roiNames = vfc.list_roiNames;

% data types we want to look at
list_dtNames  = vfc.list_dtNames;

% names of the rm in each dt
list_rmNames  = vfc.list_rmNames;

%% define things
numRois = length(list_roiNames);
numRms = length(list_dtNames);
numSubs = length(list_subInds);

%% get the rmroi cell
if isempty(rmroiCell)
    % Calculate rmroiCell
    rmroiCell = ff_rmroiCell(cr,list_subInds, list_roiNames, list_dtNames, list_rmNames);
else
    % Check if the size is ok: 
    if ~isequal(size(rmroiCell),[length(list_subInds),length(list_roiNames),length(list_dtNames)])
        error('Size of rmroiCell is not the same as its components')
    end
end

%% make averaged coverage plot for each roi and rm model
%% loop over rois
for jj = 1:numRois
    
    
    %% this roi
    roiName = list_roiNames{jj};
    
    % loop over dts
    % GLU edit 
    for kk = 1:numRms
        % figure; 
        
        % name of this dt and rm
        dtName = list_dtNames{kk};
        rmName = list_rmNames{kk};    
        
        % counter for subjects with valid roi definitions
        counter = 0;
                
        % dtName = list_dtNames{kk};
        % R = rmroi(kk,:);
        % figure; 
        mrvNewGraphWin([roiName '- ' rmName],[],visibility);
        
        RF_mean = ff_rmPlotCoverageGroup(rmroiCell(:,jj,kk), ...
                     vfc, ...
                     'flip',flip, ...
                     'visibility',visibility, ...
                     'fname', [roiName '- ' rmName]);
        
        % set user data to have RF_mean
        set(gcf, 'UserData', RF_mean);
        
        %% plot the contour if so desired
        if plotContour
            [contourMatrix, contourCoordsX, contourCoordsY] = ...
                 ff_contourMatrix_makeFromMatrix(RF_mean,vfc,contourLevel); 

            % transform so that we can plot it on the polar plot
            % and so that everything is in units of visual angle degrees
            contourX = contourCoordsX/vfc.nSamples*(2*vfc.fieldRange) - vfc.fieldRange; 
            contourY = contourCoordsY/vfc.nSamples*(2*vfc.fieldRange) - vfc.fieldRange; 
            
            % plotting
            plot(contourX, contourY, ':', 'LineWidth',2, 'Color', [0 0 0])
        end
        
        axis([-vfc.fieldRange vfc.fieldRange -vfc.fieldRange vfc.fieldRange])
        
        %% save 
        titleName = {
            ['Group Avg Coverage'], ...
            [roiName '- ' rmName], ...
            ['vfc.method: ' vfc.method]
            };
        title(titleName, 'FontWeight', 'Bold')
        if savefig
            % Save the figure
            saveas(gcf,fullfile(cr.dirs.FIGPNG, ...
                strcat([fname strrep(strjoin(titleName),' ','_') '_' vers],'.png')),'png');
            saveas(gcf,fullfile(cr.dirs.FIGSVG, ...
                strcat([fname strrep(strjoin(titleName),' ','_') '_' vers],'.svg')),'svg');
        end
                
    end
    
end


end
