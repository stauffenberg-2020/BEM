function [General, op_pts, Blade] = read_turbine_file(file)

    fid = fopen(file,'rt');
    tline = fgetl(fid);
    i=1;
    while ischar(tline)
        filebyline{i,1} = tline;
        tline = fgetl(fid);
        i=i+1;
    end
    fclose(fid);
    
    header_rows = strncmpi('#',filebyline,1);
    header_row_ids = find(header_rows);
    header_titles = filebyline(header_row_ids);
    for i=1:length(header_titles)
        header_titles(i,:) = extractAfter(header_titles(i),'# ');
        header_titles(i,:) = strtrim(header_titles(i));
    end

    for i=1:length(header_titles)
        switch header_titles{i}
            case 'GENERAL'
                General.N = str2double(extractBefore((filebyline(header_row_ids(i)+1)),';'));
                General.rho = str2double(extractBefore((filebyline(header_row_ids(i)+2)),';'));
                General.induction = str2double(extractBefore((filebyline(header_row_ids(i)+3)),';'));
                General.tip_loss = str2double(extractBefore((filebyline(header_row_ids(i)+4)),';'));
                General.highCT = str2double(extractBefore((filebyline(header_row_ids(i)+5)),';'));
            case 'OPERATIONAL_SET-POINTS'
                headers = split_string(filebyline{header_row_ids(i)+1},' ');
                data = zeros(header_row_ids(i+1)-header_row_ids(i)-3,length(headers));
                for j=1:header_row_ids(i+1)-header_row_ids(i)-3
                    tline = split_string(filebyline{header_row_ids(i)+j-1+2},' ');
                    data(j,1:length(headers)) = str2double(tline(1:length(headers)));
                end
                for k=1:length(headers)
                    op_pts.(headers{k})=data(:,k);
                end
            case 'BLADE_DETAILS'
                headers = split_string(filebyline{header_row_ids(i)+1},' ');
                data_bld = zeros(header_row_ids(i+1)-header_row_ids(i)-3,length(headers));
                for j=1:header_row_ids(i+1)-header_row_ids(i)-3
                    tline = split_string(filebyline{header_row_ids(i)+j-1+2},' ');
                    data_bld(j,1:length(headers)) = str2double(tline(1:length(headers)));
                end
                for k=1:length(headers)
                    Blade.(headers{k})=data_bld(:,k);
                end
            case 'AEROFOIL_FILE'
                aerofoil_file = strtrim(filebyline{header_row_ids(i)+1});
                [Blade.pro_t_C, Blade.pro_AoA, Blade.pro_cL, Blade.pro_cD] = read_pc_file(aerofoil_file);
        end
    end
end

