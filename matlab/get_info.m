function [beta,g,g_err]=get_info(x)
    f_address=sprintf('/Users/Launch/Coding Projects/Research/Matlab_code/matlab/L%g_c0.30_ssc.dat.new3',x)
    fid = fopen(f_address);
    dat = textscan(fid, '%f %f %f %f %f %f');
    beta = dat{4};
    g = dat{5};
    g_err = dat{6};
    fclose(fid);
end