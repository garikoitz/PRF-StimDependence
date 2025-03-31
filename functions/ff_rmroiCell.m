function [rmroiCell, selindexCell] = ff_rmroiCell(cr, ...
                                                  list_subInds, ...
                                                  list_roiNames, ...
                                                  list_dtNames, ...
                                                  list_rmNames, ...
                                                  varargin)
% rmroiCell = ff_rmroiCell(list_subInds, list_roiNames, list_dtNames, list_rmNames, varargin)
% 
% If not specified, list_path will be list_sessionRet by default
% Otherwise, specificy the key value pair, 
% e.g.:   ff_rmroiCell( ...., 'list_path', list_sessionSizeRet)
% list_sessionSize
%
%
% Makes a cell of rmrois.
% rmroiCell is a i x j x k cell, where the i corresponds to subject, j
% corresponds to ROI, and k corresponds to ret model
% Often one will call ff_rmroiGetSameVoxels(rmroiCell) afterwards
% Functionalized for code readability
% ---------------------------------------
% TSeries related calculations are slow
% the findPeaks function cannot be done on the matrix so is slow for loops
% only do if necessary ...
calcPeaks = false;
calcTSeries = false; 

%% Deal with optional input arguments
% bookKeeping; % default values

p = inputParser; 
p.addRequired('cr'            , @isstruct);
p.addRequired('list_subInds'  , @isfloat);
p.addRequired('list_roiNames' , @iscell);
p.addRequired('list_dtNames'  , @iscell);
p.addRequired('list_rmNames'  , @iscell);
p.addParameter('list_path', cr.bk.list_sessionRet);
p.addParameter('latest_fFit', false, @islogical)
p.addParameter('checkYear', '2022', @ischar)
p.addParameter('data_type', 'rmroiCell', @ischar)
% Parse it
p.parse(cr, list_subInds, list_roiNames, list_dtNames, list_rmNames, varargin{:});
% Assign it
list_path         = p.Results.list_path; 
latest_fFit       = p.Results.latest_fFit; 
checkYear         = p.Results.checkYear; 
data_type         = p.Results.data_type; 

% Validate
rm_roi = false; time_s = false;
if strcmp(data_type, 'rmroiCell')
    rm_roi = true;
elseif strcmp(data_type, 'TS')
    time_s = true;
else
    error([data_type 'not known'])
end

%% Define things
numSubs = length(list_subInds);
numRois = length(list_roiNames);
numRms  = length(list_rmNames);

rmroiCell    = cell(numSubs, numRois, numRms);
selindexCell = cell(numSubs, numRois);

%% Do things

display('Making the rmroi cell  -----------')

for ii = 1:numSubs
        
   % clear global
   mrvCleanWorkspace;
    
   subInd =  list_subInds(ii);
   subInitials = cr.bk.list_sub{subInd};
   dirVista = list_path{subInd};
   dirAnatomy = cr.bk.list_anatomy{subInd};
   chdir(dirVista);
   setpref('VISTA', 'fileFormat', 'nifti')
   load(fullfile(dirVista,'mrSESSION.mat'))
   t1_file_path = fullfile(dirAnatomy,'t1.nii.gz');
   setVAnatomyPath(t1_file_path);

   vw = initHiddenGray;

   for kk = 1:numRms

       dtName = list_dtNames{kk};
       rmName = list_rmNames{kk};
       if latest_fFit
           rmBasePath = fullfile(dirVista,'Gray',dtName);
           allFiles = dir(fullfile(rmBasePath,['retModel-' dtName '-css-fFit*']));
           [~,idx]  = sort([allFiles.datenum],'descend');
           if length(idx) > 0
                latest   = allFiles(idx(1));
                % Check it is 2022 just in case
                assert(contains(latest.date,checkYear))
                rmPath = fullfile(latest.folder,latest.name);
                fprintf('\n Will read %s, created: %s \n\n',rmPath,latest.date)
           else
               rmPath = '';
           end
       else
           rmPath = fullfile(dirVista,'Gray',dtName, rmName);
       end
       rmExists = exist(rmPath,'file');
       tsPath = fullfile(dirVista,'Gray',dtName,'TSeries','Scan1','tSeries1.mat');
       tsExists = isfile(fullfile(dirVista,'Gray',dtName,'TSeries','Scan1','tSeries1.mat'));
       
       vw = viewSet(vw, 'curdt', dtName); 
       
       % load the ret model if it exists
       if rmExists
           vw = rmSelect(vw, 1, rmPath);
           vw = rmLoadDefault(vw);
       end

       % If we want TS, load them here
       if time_s && tsExists
           allTS = load(tsPath);
           allTS = allTS.tSeries;
       end

       for jj = 1:numRois

           roiName = list_roiNames{jj};
           roiPath = fullfile(dirAnatomy, 'ROIs', [roiName '.mat']);
           [vw, roiExists] = loadROI(vw, roiPath, [],[],1,0);

           if roiExists && rmExists && tsExists

               % get the rmroi params and store it
               rmroi = rmGetParamsFromROI(vw);

               if calcTSeries
                   % add the amplitude metric
                   % first get roi coordinates and time series
                   roiCoords = viewGet(vw, 'roiCoords');
                   [tSeriesCell, ~] = getTseriesOneROI(vw,roiCoords,[], 0, 0 );
                   tSeries = tSeriesCell{1};
                   clear tSeriesCell;
                   numCoords = size(roiCoords, 2);

                   % get the mean of the top 8 values
                   numPoints = 8;
                   meanMax = ff_tSeries_meanOfMaxPoints(tSeries, numPoints);
                   rmroi.meanMax = meanMax;

                   % get the mean of the peaks
                   meanPeaks = zeros(1,numCoords);
                   if calcPeaks
                       display('Calculating peak information')
                       for vv = 1:numCoords
                           pks = findpeaks(tSeries(:,vv), 'NPeaks', numPoints, ...
                               'SortStr', 'descend', ...
                               'MinPeakDistance',5);
                           meanPeaks(vv) = mean(pks);
                       end
                   end
                   rmroi.meanPeaks = meanPeaks;
               end
               if rm_roi
                   rmroiCell{ii,jj,kk} = rmroi;
               end
               if time_s
                   % Filter to the ROI
                   TS = allTS(:,rmroi.indices);
                   rmroiCell{ii,jj,kk} = TS;
               end
           else
               % if we come here, either the roi or the rm does not exist.
               % Be informative.
               if ~roiExists, display([roiName 'for ' subInitials ' does not exist. ']), end
               if ~tsExists, display(['TS file for ' subInitials ' does not exist. ']), end
               if ~rmExists, display([rmName 'for ' subInitials ' does not exist. ']), end
               % assign D{ii,jj,kk} to be the empty cell
               rmroiCell{ii,jj,kk} = [];
           end

       end

   end
end

end