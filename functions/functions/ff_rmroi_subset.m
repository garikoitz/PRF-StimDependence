function rmroiCoord = ff_rmroi_subset(rmroi, indx)
% rmroiCoord = ff_rmroi_getSingleCoord(rmroi, indx)
%
%% for a rmroi struct consisting of rm information for numVoxels, just grab
% the information for a single voxel (as indicated by ind)

rmroiCoord = rmroi; 
if isfield(rmroi, 'coords')
    if ~isempty(rmroi.coords)
        rmroiCoord.coords = rmroi.coords(:,indx);
    end
end
if isfield(rmroi, 'indices')
    rmroiCoord.indices = rmroi.indices(indx);
end
rmroiCoord.co = rmroi.co(indx); 
rmroiCoord.sigma1 = rmroi.sigma1(indx); 
rmroiCoord.sigma2 = rmroi.sigma2(indx); 
% rmroiCoord.theta = rmroi.theta(indx); 
% rmroiCoord.beta = rmroi.beta(indx, :);
rmroiCoord.x0 = rmroi.x0(indx);
rmroiCoord.y0 = rmroi.y0(indx);
rmroiCoord.sigma = rmroi.sigma(indx);
% rmroiCoord.exponent = rmroi.exponent(indx);
rmroiCoord.polar = rmroi.polar(indx);
% rmroiCoord.rawrss = rmroi.rawrss(indx);
% rmroiCoord.rss = rmroi.rss(indx);
% rmroiCoord.thetaCenters = rmroi.thetaCenters(indx);

% NASTY GLU FIX IT
% rmroiCoord.ph = rmroi.ph(indx);
rmroiCoord.ph = zeros(size(rmroi.ecc(indx)));
rmroiCoord.ecc = rmroi.ecc(indx);
% rmroiCoord.betaScale = rmroi.betaScale(indx);

% these fields are not always computed because it takes a while
if isfield(rmroi, 'betaScale')
    if ~isempty(rmroi.betaScale)
        rmroiCoord.betaScale = rmroi.betaScale(indx);
    end
end
if isfield(rmroi, 'meanMax')
    if ~isempty(rmroi.meanMax)
        rmroiCoord.meanMax = rmroi.meanMax(indx);
    end
end

if isfield(rmroi, 'meanPeaks')
    if ~isempty(rmroi.meanPeaks)
        rmroiCoord.meanPeaks = rmroi.meanPeaks(indx);
    end
end

end