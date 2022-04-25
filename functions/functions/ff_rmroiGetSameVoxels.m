function rmroiCellSameVox = ff_rmroiGetSameVoxels(rmroiCell, vfc)
%
% sameVoxRmroi = ff_rmroiGetSameVoxels(rmroiCell, vfc)
% We want to be sure that we're looking at the same voxels for each ROI,
% across all ret models.
%
% GLU 08-2020 Edit so that we have the pRFs in the right hemifield
%
% INPUTS:
% - rmroiCell (cell array of rmrois (corresponding to different stimuli, e.g.) for a subject)
% - vfc
%   vfc.cothresh
%   vfc.eccthresh
%   vfc.sigmaMajthresh
%   vfc.sigmaEffthresh
%
% OUTPUTS: 
% sameVoxRmroi: a cell of the same dimension as rmroiCell. Each element is an
% rmroi struct
%
%%
% initialize
rmroiCellSameVox = cell(size(rmroiCell)); 

% get a single rmroi to get basic info
rmroi = rmroiCell{1}; 

    if ~isempty(rmroi)

    % the number of voxels in the original roi definition
    numVoxels = length(rmroi.indices); 

    % number of ret models we are working with
    numRms = length(rmroiCell); 

    % initialize
    indxMaster = ones(1,numVoxels); 

    %% indices bookKeeping
    for kk = 1:numRms

        rmroi = rmroiCell{kk}; 

        % the indices that pass the varExp threshold
        coindx = (rmroi.co > vfc.cothresh);

        % the indices that pass the ecc threshold
        eccindx = (rmroi.ecc < vfc.eccthresh(2)) & (rmroi.ecc > vfc.eccthresh(1));

        % sigma major threshold
        sigmaMajindx = (rmroi.sigma1 < vfc.sigmaMajthresh(2)) & ...
            (rmroi.sigma1 > vfc.sigmaMajthresh(1)); 
        
        % sigma effective threshold
        sigmaEffindx = (rmroi.sigma < vfc.sigmaEffthresh(2)) & ...
            (rmroi.sigma > vfc.sigmaEffthresh(1)); 
        
        % the indices below a certain co threshold
        % useful for when we are looking for noise
        % specify 1 otherwise
        coceilindx = (rmroi.co <= vfc.cothreshceil);
        
        % The center needs to be inside this quadrant
        % Clarify here what polar is measuring
        % [P,ECC] = cart2pol(rmroi.x0,rmroi.y0);
        % rmroi.ph == P 
        % histogram(rad2deg(rmroi.polar));hold on
        % histogram(rad2deg(rmroi.ph));hold on
        % I need ph, rads going ccw, not Polar. same but cw
        
        %{
        quads    = floor(rad2deg(rmroi.ph)/90)+1;
        hemiindx = ismember(quads, vfc.quadthresh);     
        
        % the indices after looping through rmrois
        indx = coindx & eccindx & sigmaMajindx & sigmaEffindx ... 
               & coceilindx & hemiindx; 
        %}   
        indx = coindx & eccindx & sigmaMajindx & sigmaEffindx & coceilindx;
           
           

        indxMaster = indx & indxMaster; 
        
    end

        %% the same voxels
        for kk = 1:numRms
            rmroi = rmroiCell{kk}; 
            newRmroi = ff_rmroi_subset(rmroi, indxMaster); 
            rmroiCellSameVox{kk} = newRmroi;  
        end
    
    end
end