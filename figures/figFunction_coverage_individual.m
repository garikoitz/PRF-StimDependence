function [RFcov,weight, data] = figFunction_coverage_individual(cr, opt, varargin)
%Reading circuitry field of view paper (2016).
% Makes figures related to individual visual field coverage
% Includes options for contours and ellipses

% This was a script created by Rosemary Le, and part of coverageReading
% repository. For the Stimulus Dependence paper (2019) we tried to make it into
% a reproducible and reusable process. The thing is that we need to separate
% this specific project from others, we need to separate data from code. 


%% Read parameters (this is usually the modify section of RL's scripts)
%% modify here
% This content has been moved to two places:
%     cr.codeDir/defineDefaults/IndivFigDefaults.m
%     cr.dataDir/matlabFiles/defineProjectDefaults.m

% Make varargin lower case, remove white spaces...
varargin = mrvParamFormat(varargin);
% Parse
p = inputParser;
p.addRequired('cr'  , @isstruct);
p.addRequired('opt' , @isstruct);

% Parse. Assign result inside each case
p.parse(cr, opt, varargin{:});
% Read here only the generic ones
% opt = p.Results.opt;


%%

    
    dirVista    = opt.list_path{opt.list_subInds}; 
    dirAnatomy  = cr.bk.list_anatomy{opt.list_subInds};
    subInitials = cr.bk.list_sub{opt.list_subInds};
    chdir(dirVista); 
    vw          = initHiddenGray; 
    
    for kk = 1:length(opt.list_dtNames)
        
        % load the ret model
        dtName     = opt.list_dtNames{kk}; 
        rmName     = opt.list_rmNames{kk};
        rmPath     = fullfile(dirVista, 'Gray', dtName, rmName); 
        rmExists   = exist(rmPath, 'file');
        rmDescript = opt.list_rmDescripts{kk};
        
        for jj = 1:length(opt.list_roiNames)
            
            % load the roi
            roiName = opt.list_roiNames{jj}; 
            roiPath = fullfile(dirAnatomy, 'ROIs', [roiName '.mat']);
            [vw, roiExists] = loadROI(vw, roiPath, [], [], 1, 0);
            
            % if roi and ret model exists ...
            if rmExists && roiExists
                
                % load the ret model
                vw = rmSelect(vw, 1, rmPath);
                vw = rmLoadDefault(vw); 
                
                % get the rmroi
                rmroi = rmGetParamsFromROI(vw); 
                
                % plot!
                [RFcov, ~, ~, weight, data] = rmPlotCoveragefromROImatfile(rmroi,opt.vfc);
                
                % Info for plotting purposes
                contourLevel = opt.vfc.ellipseLevel; 
                roiNameDescript = ff_stringRemove(roiName, '_rl');                     
                infoString = [roiNameDescript '. ' rmDescript '. Contour ' num2str(contourLevel) '. ' subInitials];
                    
                % title
                titleName = {
                    opt.titleDescript;
                    infoString; 
                    ['vfc.method: ' opt.vfc.method]
                    mfilename;
                    };
                title(titleName, 'FontWeight', 'Bold', 'FontSize',13)

                set(gcf, 'Color', 'w');
                
                % save
                % ff_dropboxSave; 

                %% Plot just the ellipse                         
                %% Plot the ellipse and the contourLevel
                %% Plot the ellipse over the max profile coverage
     
            end % if roi and ret model exists
            
        end % loop over rois
        
    end % loop over ret models
    
end

