function [RF_mean, RF_individuals,thisEmpties] = ff_rmPlotCoverageGroup(M, vfc, varargin)
%% function to plot the group average visual field coverage
% [RF_mean, RF_individuals] = ff_rmPlotCoverageGroup(M, vfc)
% [RF_mean, RF_individuals] = ff_rmPlotCoverageGroup(M, vfc, 'flip', false)
% INPUTS
% M:    M should be a cell vector of length numSubs  where M{ii} is the rm struct for a
        % particular roi for the iith subject
% vfc   visual field coverage information thresholds. TODO: provide more detail here. 
%
% OUTPUTS
% RF_mean: a 128 x 128 matrix
% RF_individuals: a 128 x 128 x numSubs matrix

% NOTE: there are 2 places where we do y flipping: lines 35 andn 57

% RL 03/15
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of subjects to average over
numSubs     = length(M); 
numToAvg    = numSubs; 

% % store each subject's rf (prior to averaging)
% RF      = zeros(128,128, numSubs); 
% 
% % initialize the mean visual field coverage
% RF_mean = zeros(128,128); 
RF = []; 
RF_mean = [];

% initalize number of valid subjects
counter = 0; 

% initialize vector: centers across all subjects
centersx0 = []; 
centersy0 = []; 


%% input parser -- flipping
par = inputParser; 

par.addRequired('M'        , @iscell);
par.addRequired('vfc'      , @isstruct);

addOptional(par, 'flip'      , true, @islogical); 
addOptional(par, 'fname'     , ''  , @ischar); 
addOptional(par, 'visibility', 'on', @ischar); 
addOptional(par, 'density', false, @islogical);
addOptional(par, 'size_factor', 1, @isnumeric);

parse(par, M, vfc, varargin{:});
flip       = par.Results.flip; 
fname      = par.Results.fname; 
visibility = par.Results.visibility; 
density    = par.Results.density; 
size_factor= par.Results.size_factor;

%% get the rf information for each subject
thisEmpties = [];
for ii = 1:numSubs
   if ~isempty(M{ii})
        [rf, figHandle, all_models, weight, data] = ...
            rmPlotCoveragefromROImatfile(M{ii}, vfc, size_factor);
        
        % y flippage!
        if flip
            rf = flipud(rf);
            if ~isempty(data)
                data.suby0 = -data.suby0;
            end
        end

                
        % check that rf isn't all nans. (this happens when we don't have a ret model or roi defined?)
        coverageHasNans = (sum(sum(~isnan(rf))) == 0);
        
        % Update 04/2017. If the RF is all nans it should not be averaged.
        % But if it is all zeros, it should be ...
        
        % also check that there voxels pass the threshold!!! (otherwise the coverage will be all 0s)
        % if no voxels pass threshold, <data> from rmPlotCoveragefromROImatfile is empty
        
        % Brian and I decided that if no voxel passes, then it should not
        % average, it is going to change the N, that's it. 
        
        noVoxelsPassThreshold = isempty(data);               
        
        if ~coverageHasNans && ~noVoxelsPassThreshold 
        % if ~coverageHasNans
            
            counter = counter + 1;


            if density
                warning('We using avg and density true')
                rf = rf ./ length(data.suby0);
                rf = rf ./ max(rf(:));
            end


            % store it in RF
            RF(:,:,counter) = rf; 
            if ~noVoxelsPassThreshold
                close
            end
        else
            ii
            thisEmpties = [thisEmpties, ii];
        end
        
        % grab location of rf centers
        if ~isempty(data)
            centersx0 = [centersx0, data.subx0]; 
            centersy0 = [centersy0, data.suby0];
        end
        
   end
end

if ~isempty(RF)

    %% and average them together ...

    RF_mean = mean(RF,3); 

    % flip about the x axis
    RF_mean = flipud(RF_mean);
    centersy0 = -centersy0;

    % plot the average coveraged
    % edit GLU
    % figure();
    mrvNewGraphWin(fname,[],visibility);
    
    % to make the black outer circle thing
    c = makecircle(vfc.nSamples);  
    % to make the polar angle plot
    inc = linspace(-vfc.fieldRange,vfc.fieldRange, vfc.nSamples);


    RF_mean = RF_mean.*c; % make black outer circle
    imagesc(inc,inc',RF_mean); 

    % add polar grid on top
    p.ringTicks = (1:3)/3*vfc.fieldRange;
    p.color = 'k';
    
    % tick labels
    if ~isfield(vfc, 'tickLabel')
        p.tickLabel = true; 
    else
        p.tickLabel = vfc.tickLabel;
    end
    
    polarPlot([], p);

    title(['Group Average.'], 'FontSize', 14);
    axis off

    %% colorbar things

    % when averaging, values don't necessarily take the full range of 0 - 1
    % but we want the average plots to look the same, so we have the colormap
    % range from 0-1
    c = colorbar;
    colormap(vfc.cmap);
    caxis([0 1]);

    %% add the centers if requested
    if vfc.addCenters
        plot(centersx0, centersy0,'.','Color',[.5 .5 .5], 'MarkerSize', 4);
    end
 
    %% group summary contour 
    if vfc.contourPlot

        % get the x and y points of the specified contour level
        [~, contourCoordsX, contourCoordsY] = ...
        ff_contourMatrix_makeFromMatrix(RF_mean,vfc,vfc.contourLevel);

        % transform so that we can plot it on the polar plot
        contourX = contourCoordsX/vfc.nSamples*(2*vfc.fieldRange) - vfc.fieldRange; 
        contourY = contourCoordsY/vfc.nSamples*(2*vfc.fieldRange) - vfc.fieldRange;
        
        plot(contourX, contourY, '--', 'Color', vfc.contourColor, 'MarkerSize', 2, 'LineWidth',2);
        
    end
    
end

RF_individuals = RF; 

end