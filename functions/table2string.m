function tstr = table2string(T, valid_sessions)
%table2string converts table to a printable from in a string
%   Obtains subjects and sessions summary

%{
  sub = categorical({'sub-07','sub-08','sub-08','sub-09','sub-10'}');
  ses = categorical({'ses-09','ses-07','ses-08','ses-09','ses-09'}');
  valid_sessions = logical([1,0,1,0,1]');
  T = table(sub, ses, valid_sessions);
%}
    T.valid_sessions = valid_sessions;

    tstr = [];
    subs = unique(T.sub);
    
    for ns = 1:length(subs)
        
        sess = unique(T.ses(T.sub==subs(ns)));
        ses_str = '(';
        for nn=1:length(sess)
            is_valid = T.valid_sessions(T.sub==subs(ns) & T.ses==sess(nn));
            ses = strip(strrep(char(sess(nn)), 'ses-',''),'left','0');
            if is_valid
                ses = ['\bf' ses '\rm'];
            end
            ses_str = [ses_str ses ','];
        end
        tstr = [tstr char(subs(ns)) strip(ses_str,'right',',') '); '];
    end
    tstr = strip(tstr,'right',' ');
    tstr = strip(tstr,'right',';');
end




