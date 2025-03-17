function sd_create_distr_plot(  same_vox_data, ...
                                data_names, ...
                                roi_names, ...
                                R2, ...
                                varargin)

%{
same_vox_data = CNI_R;
data_names = CNI_data_types;
roi_names = new_list_roiNames;
R2 = varExplained;
show_fig = 'off';
save_fig = true;
path_fname = path_fname;
filter_value = 5;
%}

%% Varargin
% Make varargin lower case, remove white spaces...
varargin = mrvParamFormat(varargin);
% Parse
p = inputParser;
p.addRequired('same_vox_data', @isstruct);
p.addRequired('data_names', @iscell);
p.addRequired('roi_names', @iscell);
p.addRequired('R2', @isnumeric);

p.addOptional('save_fig' , false , @islogical);
p.addOptional('path_fname' , '' , @ischar);
p.addOptional('filter_value' , 5 , @isnumeric);
p.addOptional('show_fig' , 'off' , @ischar);

% Parse. Assign result inside each case
p.parse(same_vox_data, ...
        data_names, ...
        roi_names, ...
        R2, ...
        varargin{:});
% Read here only the generic ones
show_fig  = p.Results.show_fig;
save_fig  = p.Results.save_fig;
path_fname = p.Results.path_fname;
filter_value = p.Results.filter_value;

% Input operations
DATA = same_vox_data.rmroiCellSameVox;

if size(DATA,3)~=2
    error("Data's 3rd dim needs to be 2, we are comparing two distributions.")
end
if size(data_names, 2) ~= size(DATA,3)
    error("Number or data names doesn't not correspond to the data size")
end
if size(roi_names, 1) ~= size(DATA,2)
    error("Number or ROI names doesn't not correspond to the data size")
end
fprintf("Size of data is %i subjects and %i ROIs.\n", size(DATA,1), size(DATA,2))

if save_fig && length(path_fname)==1
    error("Give a path and a name if you want to save the file")
end

% Do the plotting
[pathname,fname,fext]=fileparts(path_fname);
f = mrvNewGraph(fname,'wide',show_fig);

nrows = 2; ncols = ceil(length(roi_names)/nrows);
ii=0;
for nr=1:length(roi_names)
% for nr=7:8
    ii = ii + 1;
    subplot(nrows, ncols, ii)

    roi = roi_names{nr};
    data1=[];
    data2=[];
    data1_filt=[];
    data2_filt=[];
    N = [];
    N_filt = [];
    y_maxs = [];
    for ns=1:size(DATA,1)
        data1 = [data1, DATA{ns, nr, 1}.ecc];
        data2 = [data2, DATA{ns, nr, 2}.ecc];
        data1_filt=[data1_filt, data1(data2>filter_value)];
        data2_filt=[data2_filt, data2(data2>filter_value)];
        N = [N, length(data1)];
        N_filt = [N_filt, length(data1_filt)];
    end
    
    [x_values1, mu1, sigma1, mn1, mx1] = dr_distPlottingVals(data1, 90);
    [x_values2, mu2, sigma2, mn2, mx2] = dr_distPlottingVals(data2, 90);
    [x_values1_filt, mu1_filt, sigma1_filt, mn1_filt, mx1_filt] = dr_distPlottingVals(data1_filt, 90);
    [x_values2_filt, mu2_filt, sigma2_filt, mn2_filt, mx2_filt] = dr_distPlottingVals(data2_filt, 90);
    
    [indiv_pdf1, x_values1] = ksdensity(data1);
    [indiv_pdf2, x_values2] = ksdensity(data2);
    [indiv_pdf1_filt, x_values1_filt] = ksdensity(data1_filt);
    [indiv_pdf2_filt, x_values2_filt] = ksdensity(data2_filt);
    mySupTitle = 'Density distribution of eccentricity values.';
    at1  = plot(x_values1, indiv_pdf1, 'Color', 'r', ...
                                                    'LineWidth',1.5, 'LineStyle',':');
    hold on;
    at2  = plot(x_values2, indiv_pdf2, 'Color', 'b', ...
                                                    'LineWidth',1.5, 'LineStyle',':');
    at1_filt  = plot(x_values1_filt, indiv_pdf1_filt, 'Color', 'r', ...
                                                    'LineWidth',2, 'LineStyle','-');
    at2_filt  = plot(x_values2_filt, indiv_pdf2_filt, 'Color', 'b', ...
                                                    'LineWidth',2, 'LineStyle','-');
    
    % Calculate ymax
    y_max = max([max(indiv_pdf1), max(indiv_pdf2), max(indiv_pdf1_filt), max(indiv_pdf2_filt)]);
    xlim([0, 15]); ylim([0, y_max]);
    % xticks([0:0.2:.8]);yticks([0:5:20]);
    % title(tn)
    a = [at1_filt;at2_filt];  % Concatenate plots to be used with legend
    % l = [l; strcat(strrep(cats{sc},'_','\_'),'TEST'); ...
    %         strcat(strrep(cats{sc},'_','\_'),'RETEST')];
    % Calculate Cohen's d
    d = (mu2 - mu1) / sqrt((sigma2^2 + sigma1^2) / 2);
    d_filt = (mu2_filt - mu1_filt) / sqrt((sigma2_filt^2 + sigma1_filt^2) / 2);

    % Add vertical lines to indicate mean and Cohen's d
    text(.25, y_max - 0.02, sprintf("%s, R^2>%.1g", strrep(strrep(roi,"WangAtlas_",""),"_","\_"), R2))
    xline(mu1, 'Color', 'r', 'LineStyle', ':','LineWidth',.5);
    xline(mu2, 'Color', 'b', 'LineStyle', ':','LineWidth',.5);
    text(mu1, (y_max - 0.05), ...
         sprintf("d':%1.2g, N:%4.1f(%4.1f)", d, mean(N),std(N)),...
            'LineStyle','-')
    xline(mu1_filt, 'Color', 'r', 'LineStyle', '-');
    xline(mu2_filt, 'Color', 'b', 'LineStyle', '-');
    text(mu1_filt, .06, ...
        sprintf("d'_{>%i}: %1.2g, N:%4.1f(%4.1f))", ...
                filter_value, ...
                d_filt, ...
                mean(N_filt), ...
                std(N_filt)), ...
        'LineStyle','-')
    legend(a, data_names)
    xlabel('\delta eccentricy')
end
% Save figure
if save_fig
    fprintf("Saving file as %s.\n", path_fname)
    if isfile(path_fname)
        delete(path_fname)
    end
    if ~isdir(pathname)
        mkdir(pathname)
    end
    saveas(f, path_fname);
end

end
