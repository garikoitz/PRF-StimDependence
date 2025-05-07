function [R,C_data,cr, icc_table]=sd_filter_same_voxel_pairs(cr,...
                                         rmroiCell, ...
                                         list_subInds,...
                                         list_roiNames,...
                                         vfc, ...
                                         varargin)
%% Varargin
% Make varargin lower case, remove white spaces...
varargin = mrvParamFormat(varargin);
% Parse
p = inputParser;
p.addRequired('cr'            , @isstruct);
p.addRequired('rmroicell'     , @iscell);
p.addRequired('list_subInds');
p.addRequired('list_roiNames' , @iscell);
p.addRequired('vfc'           , @isstruct);
p.addOptional('show_summary'  , false , @islogical);
p.addOptional('calculate_icc' , false , @islogical);
p.addOptional('rmNames'       , {'NA', 'NA'}, @iscellstr);
p.addOptional('icc_meas'      , 'ecc', @ischar);


subistable = false;
if istable(list_subInds); subistable = true; end

% Parse. Assign result inside each case
p.parse(cr, rmroiCell, list_subInds, list_roiNames, vfc, varargin{:});
% Read here only the generic ones
show_summary  = p.Results.show_summary;
calculate_icc = p.Results.calculate_icc;
rmNames       = p.Results.rmNames;
icc_meas      = p.Results.icc_meas;

% INITIALIZE SOME THINGS
numRois = length(list_roiNames);
if subistable
    numSubs = height(list_subInds);
else  
    numSubs = length(list_subInds);
end
% if calculate_icc is true, see if a table of subs was passed and only one
% sub was passed. 
icc_table = table();
if calculate_icc
    if ~subistable
        warning(['calculate_icc=true will be ignored, as list_subInds is not '...
                'a table. This will be only used to calculate the ICC for one '...
                'subject across different sessions'])
    end
    if length(unique(list_subInds.sub)) > 1
        warning(['calculate_icc=true will be ignored, as list_subInds is a table ' ...
                 'with more than one subject. This will be only used to calculate ' ...
                 'the ICC for one subject across different sessions'])
    end
    disp('Conditions to calculate per-subject intra-session ICC seem ok, doing')
    roi_data_icc = NaN(1, numRois);
    for nm = 1:2
        if ~isempty(rmroiCell{1,1,nm})
            for nr = 1:numRois
                M = NaN(numSubs, size(rmroiCell{1,nr,nm}.(icc_meas),2));
                for ii=1:height(list_subInds)
                    if ~isempty(rmroiCell{ii,nr,nm})
                        M(ii,:) = rmroiCell{ii,nr,nm}.(icc_meas);
                    end
                end
                % Invert it so that sessions are observations
                MM = M';
                % Remove columns with NaN (not acquired or any other problem)
                MM(:, any(isnan(MM), 1)) = [];
                % If there is only one sub, ignore, no ICC
                if size(M,2) > 1
                    roi_data_icc(nr) = compute_icc31(MM);
                else
                    roi_data_icc(nr) = 0;
                end
            end
            % Create the table
            CT = [{unique(list_subInds.sub)}, {categorical(rmNames(nm))}, ...
                {categorical(string(icc_meas))}, num2cell(roi_data_icc)];
            T = cell2table(CT, 'VariableNames', ['sub', 'type', 'meas', list_roiNames]);
            icc_table = [icc_table; T];
        end
   end

end
% cell for linearizing the data (a vector for each ROI)
L_data  = cell(1, numRois);
X_rm1   = cell(1, numRois);
Y_rm1   = cell(1, numRois);
X_rm2   = cell(1, numRois);
Y_rm2   = cell(1, numRois);
Ecc_rm1 = cell(1, numRois);
Ecc_rm2 = cell(1, numRois);

% In comparing ret models, the collection of voxels may not be the same
% because of the thresholding. In this cell we redefine the rmroi
rmroiCellSameVox = cell(size(rmroiCell));
for jj = 1:numRois
    mean_subjects = NaN([numSubs,1]);
    std_subjects = NaN([numSubs,1]);
    max_subjects = NaN([numSubs,1]);
    pct_subjects = NaN([numSubs,1]);
    N_subjects = NaN([numSubs,1]);
    for ii = 1:numSubs
        % get identical voxels for each subject's roi over all ret models
        D = rmroiCell(ii,jj,:);
        if (~isempty(D{1,1,1}) && ~isempty(D{1,1,2})) 
            % GLU EDIT function: remove voxels from the oppossite hemifield
            rmroiCellSameVox(ii,jj,:) = ff_rmroiGetSameVoxels(D, vfc);
            % Calculate eccentricity differences
            ecc_diffs = rmroiCellSameVox{ii,jj,2}.ecc - rmroiCellSameVox{ii,jj,1}.ecc;
            % Calculate N
            N = length(ecc_diffs);
            % Calculate max difference and the index
            [max_val,max_ind]=max(ecc_diffs);
            mean_diff = mean(ecc_diffs);
            std_diff = std(ecc_diffs);
            % Instead of the max value, calculate a percentile value
            pctilval = 75;
            prct_diff = prctile(ecc_diffs, pctilval);
            % obtain the index value for plotting
            [~, prctile_ind] = min(abs(ecc_diffs - prct_diff));

            voxel_indices = rmroiCellSameVox{ii,jj,1}.indices;
            if show_summary
                fprintf("ROI:%s, Sub:%i, N: %i, ecc_diff_prctile-%s: %.2g (prctl_index: %i), ecc_diff_max: %.2g (max_index:%i), MEAN: %.2g, STD: %.2g\n", ...
                    list_roiNames{jj}, ...
                    ii, ...
                    N, ...
                    num2str(pctilval), ...
                    prct_diff, ...
                    voxel_indices(prctile_ind), ...
                    max_val, ...
                    voxel_indices(max_ind), ...
                    mean_diff, ...
                    std_diff);
            end
            % Add to vector to capture all subjects values
            mean_subjects(ii) = mean_diff;
            std_subjects(ii) = std_diff;
            if isempty(max_val);  max_subjects(ii) = NaN;
            else; max_subjects(ii) = max_val; end
            pct_subjects(ii) = prct_diff;
            N_subjects(ii) = N;
        else
            disp('')
        end
    end
    if show_summary
        disp('--')
        fprintf("ROI: %s, N: %.4g (%.4g), mean_prctile-%s: %.2g (%.2g), mean_max: %.2g (%.2g), mean_MEAN: %.2g (%.2g), mean_STD: %.2g (%.2g)\n", ...
                list_roiNames{jj}, ...
                mean(N_subjects), ...
                std(N_subjects), ...
                num2str(pctilval), ...
                mean(pct_subjects), ...
                std(pct_subjects), ...
                mean(max_subjects), ...
                std(max_subjects), ...
                mean(mean_subjects), ...
                std(mean_subjects), ...
                mean(std_subjects), ...
                std(std_subjects));
        disp('--')
        disp('--')
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
    cdata = ff_colormapForValues(ldata, vfc.cmapValues, vfc.cmapRange);    
    C_data{jj} = cdata;  
end


% Prepare the output
R = struct();
R.rmroiCellSameVox = rmroiCellSameVox;
R.L_data = L_data;
R.X_rm1 = X_rm1;
R.Y_rm1 = Y_rm1;
R.X_rm2 = X_rm2;
R.Y_rm2 = Y_rm2;
R.Ecc_rm1 = Ecc_rm1;
R.Ecc_rm2 = Ecc_rm2;
R.Ph_rm1 = Ph_rm1; 
R.Ph_rm2 = Ph_rm2;
R.Sm_rm1 = Sm_rm1; 
R.Sm_rm2 = Sm_rm2;

end
