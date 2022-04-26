function [s,d,m] = myCopy(src,dstf,subind,cp)
%myCopy Copy file and create dir if not there
    dstd = fileparts(dstf);
    if isfile(src)
        [s,d,m] = mkdir(dstd);
        if s 
            if cp
                st = system(['cp ' src ' ' dstf]);
                s = ~st;
            else
                [s,d,m] = copyfile(src,dstf,'f'); 
            end
            if ~s; warning('Did not work for %i, error copying to %s: %s',subind,dst,d); end
        else
            warning('Did not work for %i,could not create dir %s',subind,dstd);
        end
    else
        warning('Did not work for %i, no src file %s',subind,src);
    end
    assert(strcmp(Simulink.getFileChecksum(src),Simulink.getFileChecksum(dstf)))
end

