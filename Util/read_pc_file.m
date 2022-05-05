function [thickness, AoA, cL, cD] = read_pc_file(pc)
    fid = fopen(pc,'r');
    fid1 = fopen(pc,'r');
    
    nT = textscan(fid, '%f', 1, 'HeaderLines', 1, 'CollectOutput', 1);
    nT = nT{1,1}; % No. of tables
    entries = zeros(nT,1);
    
    j=1;
    for i=1:nT
        tline = textscan(fid, '%f%f%f', 1, 'HeaderLines', j, 'CollectOutput', 1);
        entries(i,1) = tline{1,1}(2);
        thickness(i,1) = tline{1,1}(3);
        if i==1
            section(i,1) = textscan(fid1, '%f%f%f%f%f%f', entries(i,1), 'HeaderLines', 3, 'CollectOutput', 1);
        else
            section(i,1) = textscan(fid1, '%f%f%f%f%f%f', entries(i,1), 'HeaderLines', 2, 'CollectOutput', 1);
        end
        j=entries(i,1)+1;
    end
    fclose(fid);
    fclose(fid1);
    % Linear interpolation to make sure all the data is in same size
    AoA = (-179:1:180)';
    for i=1:nT
        cL(:,i) = interp1(section{i,1}(:,1),section{i,1}(:,2),AoA);
        cD(:,i) = interp1(section{i,1}(:,1),section{i,1}(:,3),AoA);
    end
end