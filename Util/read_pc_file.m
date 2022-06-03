function [thickness, AoA, cL, cD] = read_pc_file(file)
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
    
    nT = real(str2double(tmp{2}));
    entries = zeros(nT,1);
    thickness = zeros(nT,1);
    
    AoA = (-179:1:180)';
    cL = zeros(length(AoA),length(thickness));
    cD = zeros(length(AoA),length(thickness));
    
    k=3; % Data tables start from line 3 in HAWC2 pc file
    for i=1:nT
        tline = split_string(tmp{k,1},' ');
        entries(i,1) = real(str2double(tline{1,2}));
        thickness(i,1) = real(str2double(tline{1,3}));
        section = zeros(entries(i,1),4);
        for j=1:entries(i)
            tline = split_string(strtrim(tmp{j+k,1}),' ');
            for n=1:4
                section(j,n) = real(str2double(tline{n}));
            end
        end
        % Linear interpolation to make sure all the data is in same size
        cL(:,i) = interp1(section(:,1),section(:,2),AoA);
        cD(:,i) = interp1(section(:,1),section(:,3),AoA);
        k=k+entries(i,1)+1;
    end
end