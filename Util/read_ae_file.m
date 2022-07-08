function [r, C, t_C, nT] = read_ae_file(file)
%     coder.inline('never');
    fid = fopen(file,'r');
    
    i=1;
    while ~feof(fid)
        tline = fgetl(fid);
        i=i+1;
    end
    n_lines =i-1;
    fclose(fid);
    
    fid = fopen(file,'r');
    tmp = cell(n_lines,1);

    for i=1:n_lines
        tmp{i,1} = fgetl(fid);
    end
    fclose(fid);
    
    temp = split_string(tmp{2});
    nT = real(str2double(temp{2}));
    entries = zeros(nT,1);
    thickness = zeros(nT,1);
    
    
    r = zeros(nT-1,1);
    C = zeros(nT-1,1);
    t_C = zeros(nT-1,1);
    set = zeros(nT-1,1); % Not used as of now
    
    for i=1:nT-1
        tline = split_string(tmp{3+i,1},' ');
        r(i,1) = real(str2double(tline{1,1}));
        C(i,1) = real(str2double(tline{1,2}));
        t_C(i,1) = real(str2double(tline{1,3}));
        set(i,1) = real(str2double(tline{1,4}));
    end
end