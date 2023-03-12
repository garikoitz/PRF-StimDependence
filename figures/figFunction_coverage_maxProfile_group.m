function [RF_Mean_Cells, RF_Indiv_Cells, empties] = figFunction_coverage_maxProfile_group(cr, subinds, varargin)

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

p.addOptional('flip'       , true,  @islogical);
p.addOptional('bootcontour', false,  @islogical);
p.addOptional('rmroicell' , {}, @iscell);
p.addOptional('list_roiNames' , {}, @iscell);
p.addOptional('list_dtNames' , {}, @iscell);
p.addOptional('list_rmNames' , {}, @iscell);
p.addOptional('sizedegs' , 0, @isnumeric);
p.addOptional('numboots' , 3, @isnumeric);
p.addOptional('vers','v01', @ischar);
p.addOptional('fname',''  , @ischar);
p.addOptional('invisible',false  , @islogical);
p.addOptional('minvarexp' , 0, @isnumeric);
p.addOptional('method' , 'max', @ischar);
p.addOptional('density', false  , @islogical);

% Parse. Assign result inside each case
p.parse(cr, subinds, varargin{:});
% Read here only the generic ones
flip          = p.Results.flip;
bootcontour   = p.Results.bootcontour;
rmroiCell     = p.Results.rmroicell;
vers          = p.Results.vers;
fname         = p.Results.fname;
invisible     = p.Results.invisible;
list_roiNames = p.Results.list_roiNames;
list_dtNames  = p.Results.list_dtNames;
list_rmNames  = p.Results.list_rmNames;
sizedegs      = p.Results.sizedegs;
minvarexp     = p.Results.minvarexp;
numboots      = p.Results.numboots;
method        = p.Results.method;
density       = p.Results.density;

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

if minvarexp==0
    minvarexp = vfc.cothresh; 
end
vfc.cothresh = minvarexp; 


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
if isempty(list_roiNames) 
    list_roiNames = vfc.list_roiNames;
end

% data types we want to look at
if isempty(list_dtNames)
    list_dtNames  = vfc.list_dtNames;
end

% names of the rm in each dt
if isempty(list_rmNames)
    list_rmNames  = vfc.list_rmNames;
end

% sizedegs if not passed
if sizedegs==0
    sizedegs = vfc.fieldRange;
else
    vfc.fieldRange = sizedegs;
end



%% define things
numRois = length(list_roiNames);
numRms  = length(list_dtNames);
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
RF_Mean_Cells  = cell(numRois, numRms);
RF_Indiv_Cells = cell(numRois, numRms);
empties        = cell(numRois, numRms);
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
        
        fprintf('\n\n%s >> %s >> %s\n',roiName,dtName, rmName)
        
        % counter for subjects with valid roi definitions
        counter = 0;
                
        % dtName = list_dtNames{kk};
        % R = rmroi(kk,:);
        % figure; 
        mrvNewGraphWin([roiName '- ' rmName],[],visibility);
        
        
        if bootcontour && kk==1
            numsubs    = size(rmroiCell,1);
            RFMEAN     = nan(128,128,numboots);
            for bb=1:numboots
                randReplacement = datasample(1:numsubs,numsubs);
                [RF_mean, RF_individuals, thisEmpties] = ff_rmPlotCoverageGroup(rmroiCell(:,jj,kk), ...
                    vfc, ...
                    'flip',flip, ...
                    'visibility',visibility, ...
                    'fname', [roiName '- ' rmName]);
                RFMEAN(:,:,bb) = RF_mean;
            end
        else
            [RF_mean, RF_individuals, thisEmpties] = ff_rmPlotCoverageGroup(rmroiCell(:,jj,kk), ...
                vfc, ...
                'flip',flip, ...
                'visibility',visibility, ...
                'fname', [roiName '- ' rmName]);
        end
        RF_Mean_Cells{jj,kk}  = RF_mean;
        RF_Indiv_Cells{jj,kk} = RF_individuals;
        empties{jj,kk}        = thisEmpties;
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
            
            if kk==1
                plot(contourX, contourY, ':', 'LineWidth',2, 'Color', [0 0 0])
                prevContX = contourX;
                prevContY = contourY;
            else
                plot(contourX, contourY, '--', 'LineWidth',2, 'Color', [0.5 0.5 0.5])
                if bootcontour
                    for bb=1:numboots
                        [~, contourCoordsX, contourCoordsY] = ...
                            ff_contourMatrix_makeFromMatrix(RFMEAN(:,:,bb),vfc,contourLevel);
                        contX = contourCoordsX/vfc.nSamples*(2*vfc.fieldRange) - vfc.fieldRange;
                        contY = contourCoordsY/vfc.nSamples*(2*vfc.fieldRange) - vfc.fieldRange;
                        plot(contX, contY, ':', 'LineWidth',2, 'Color', [0 0 0]); hold on
                    end
                else
                    plot(prevContX, prevContY, ':', 'LineWidth',2, 'Color', [0 0 0])
                end
            end
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
