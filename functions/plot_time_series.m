function pp = plot_time_series(pp, path_fname)

    %{
    fname = ['heb_sub01_hv4_left_ind-' pp.ind];
    fpath = '~/toolboxes/PRF-StimDependence/DATA/figures/png';
    path_fname = fullfile(fpath, [fname '.png']), 'png'
    plot_time_series(pp, path_fname)

    %}
    f = figure('visible','off');
    % f = figure('visible','on');
    
    % Extract the time series and plots them. 
    pp.input_time_series = struct();
    pp.modelpred_time_series = struct();
        % Go through all the data types pased on the struct. 
        for mn=1:length(pp.what_data_types)
            % Read the time series
            data_type_name = pp.what_data_types{mn};
            fprintf("Working with data type: %s\n", data_type_name)
            ts_all = load(fullfile(pp.p2_ret_data,data_type_name, pp.generic));
            % Select the right time series with the index and add it to the struct
            pp.input_time_series.(data_type_name) = ts_all.tSeries(:, pp.ind); 
            % Obtain the estimated time series based on the fitted parameters to this time series
            M = load(fullfile(pp.p2_ret_data, data_type_name, pp.rmNames{mn}));
            M.tSeries = pp.input_time_series.(data_type_name);
            pred = pmVistaObtainPrediction_voxel(M, pp.ind);
            pp.modelpred_time_series.(data_type_name) = pred;
            %Plot it: 
            switch  data_type_name
                case {'Words_English','Words'}
                    color = 'b';
                case {'Words_Hebrew','FalseFont'}
                    color = 'r';
                case {'Checkers'}
                    color = 'g';
                otherwise
                    color = 'k';
            end

            % plot(M.tSeries,[color '-'],'LineWidth',1); 
            
            plot(pred,[color ':'],'LineWidth',1);
            hold on; 
        end
        legend(pp.what_data_types)
    saveas(f, path_fname);
end


