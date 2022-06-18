function [pitch, lambda, out] = read_Cp_Ct(file)
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
    
    temp = split_string(tmp{1});
    nR = real(str2double(temp{1}));
    nC = real(str2double(temp{2}));
    
    lambda = zeros(1,nC);
    pitch = zeros(nR,1);
    out = zeros(nR,nC);
    
    tline = split_string(tmp{2,1},' ');
    for j=1:nC
        lambda(j) = real(str2double(tline{1,j}));
    end
    
    for i=1:nR
        tline = split_string(tmp{2+i,1},' ');
        pitch(i) = real(str2double(tline{1,1}));
        for j=1:nC
            out(i,j) = real(str2double(tline{1,j+1}));
        end
    end
end