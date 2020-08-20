function [RFcov,weight, data] = figFunction_coverage_individual(cr, subind, varargin)
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
p.addRequired('subind' , @isnumeric);

% Parse. Assign result inside each case
p.parse(cr, subind, varargin{:});
% Read here only the generic ones
% opt = p.Results.opt;


% Individual data
subname         = cr.bk.list_sub{subind};
opt             = cr.subj.(subname).params.covfig;


%%
    dirVista    = cr.bk.list_sessionRet{subind};
    dirAnatomy  = cr.bk.list_anatomy{subind};
    subInitials = cr.bk.list_sub{subind};
    
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
            % GLU ADDED THIS
            % CHeck Gomez: 
            % G = load('/share/kalanit/biac2/kgs/3Danat/jesse/ROIs/lVOTRC.mat')
            % D = load('/black/localhome/glerma/TESTDATA/PRF-StimDependence/DATA/anatomy/dames/lVOTRC.mat')
            % D2= load('/share/wandell/data/anatomy/dames/ROIs/lVOTRC.mat')
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
                
                pp = split(rmPath,filesep);
                
                
                % title
                titleName = {
                    [opt.titleDescript  ' ' num2str(subind) ' ' pp{end-4}];
                    infoString; 
                    ['vfc.method: ' opt.vfc.method]
                    mfilename;
                    };
                title(titleName, 'FontWeight', 'Bold', 'FontSize',13)

                set(gcf, 'Color', 'w');
                
                % save
                ff_dropboxSave('saveto',cr.dirs.FIG); 

                %% Plot just the ellipse                         
                %% Plot the ellipse and the contourLevel
                %% Plot the ellipse over the max profile coverage
     
            end % if roi and ret model exists
            
        end % loop over rois
        
    end % loop over ret models
    
end

