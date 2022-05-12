function [thickness, AoA, cL, cD] = read_pc_file(file)
    
    fid = fopen(file,'rt');
    tline = fgetl(fid);
    i=1;
    while ischar(tline)
        tmp{i,1} = tline;
        tline = fgetl(fid);
        i=i+1;
    end
    fclose(fid);
    
    nT = str2double(tmp{2,1});
    entries = zeros(nT,1);
    thickness = zeros(nT,1);
    
    AoA = (-179:1:180)';
    cL = zeros(length(AoA),length(thickness));
    cD = zeros(length(AoA),length(thickness));
    
    k=3; % Data tables start from line 3
    for i=1:nT
        tline = split_string(tmp{k,1},' ');
        entries(i,1) = str2double(tline{1,2});
        thickness(i,1) = str2double(tline{1,3});
        section = zeros(entries(i,1),4);
        for j=1:entries(i)
            tline = split_string(strtrim(tmp{j+k,1}),' ');
            section(j,1:4) = str2double(tline(1:4));
        end
        % Linear interpolation to make sure all the data is in same size
        cL(:,i) = interp1(section(:,1),section(:,2),AoA);
        cD(:,i) = interp1(section(:,1),section(:,3),AoA);
        k=k+j+1;
    end
end