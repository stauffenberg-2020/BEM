function [AeroTwist] = read_aero_twist(file, nsec)
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
    
    match_str = 'nsec 19';
%     plus('nsec ',string(nsec));
    for i=1:n_lines
        if strncmpi(match_str,strtrim(tmp{i}),7) == 1
            id = i;
        end
    end
    
    AeroTwist = zeros(nsec-1,1); 
    for i=1:nsec-1
        tline = split_string(tmp{id+1+i,1},' ');
        AeroTwist(i,1) = real(str2double(tline{1,6}));
    end
    AeroTwist = -AeroTwist;
end