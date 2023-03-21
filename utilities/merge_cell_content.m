function merged_cell = merge_cell_content(ccell,jj,nr,ii)
    s1 = ccell{jj, nr, ii};
    s2 = ccell{jj, nr + 3, ii};
    
    if isempty(s1) || isempty(s2)
        merged_cell = cell(1);
    else
        mname = split(s1.name,'_');
        reg = mname{2};
        reg = reg(1:2);
        newname = [mname{1} '_' reg '_' mname{3}];


        sm = s1;
        fs = fields(sm);
        for nf=1:length(fs)
            ff = fs{nf};
            switch ff
                case 'name'
                    sm.(ff) = newname;
                case {'curScan','vt','session'}
                    sm.(ff) = s1.(ff);
                case {'indices','beta'}
                    sm.(ff) = [s1.(ff); s2.(ff)];
                case {'coords','co','sigma1','sigma2','theta','x0','y0',...
                      'sigma','exponent','polar','rawrss','rss','betaScale', ...
                      'thetaCenters','ph','ecc'}    
                    sm.(ff) = [s1.(ff), s2.(ff)];
                otherwise
                    warning('Not recognized')

            end
        end

        merged_cell = {sm};
    end
end
