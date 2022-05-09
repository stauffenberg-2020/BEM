function [General, op_pts, Blade] = read_turbine_file(file)

    filestr = fileread(file);
    filebyline = regexp(filestr, '\n', 'split');
    filebyline = filebyline(:);
    filebyline = convertCharsToStrings(filebyline);

    header_rows = strncmpi('#',filebyline,1);
    header_row_ids = find(header_rows);
    header_titles = filebyline(header_row_ids);
    for i=1:length(header_titles)
        header_titles(i,:) = extractAfter(header_titles(i),'# ');
        header_titles(i,:) = strtrim(header_titles(i));
    end

    for i=1:length(header_titles)
        switch header_titles(i)
            case 'GENERAL'
                General.N = str2double(extractBefore((filebyline(header_row_ids(i)+1)),';'));
                General.rho = str2double(extractBefore((filebyline(header_row_ids(i)+2)),';'));
                General.induction = str2double(extractBefore((filebyline(header_row_ids(i)+3)),';'));
                General.tip_loss = str2double(extractBefore((filebyline(header_row_ids(i)+4)),';'));
                General.highCT = str2double(extractBefore((filebyline(header_row_ids(i)+5)),';'));
            case 'OPERATIONAL_SET-POINTS'
                headers = strsplit(filebyline(header_row_ids(i)+1));
                range_start = header_row_ids(i)+2;
                range_end = header_row_ids(i+1)-2;
                opts = detectImportOptions(file);
                opts.DataLines = [range_start range_end];
                data_tmp1 = readmatrix(file,opts);
                for j=1:length(headers)-1
                    op_pts.(headers(j))=data_tmp1(:,j);
                end
            case 'BLADE_DETAILS'
                headers = strsplit(filebyline(header_row_ids(i)+1));
                range_start = header_row_ids(i)+2;
                range_end = header_row_ids(i+1)-2;
                opts = detectImportOptions(file);
                opts.DataLines = [range_start range_end];
                data_tmp2 = readmatrix(file,opts);
                for j=1:length(headers)-1
                    Blade.(headers(j))=data_tmp2(:,j);
                end
            case 'AEROFOIL_FILE'
                aerofoil_file = strtrim(filebyline(header_row_ids(i)+1));
                [Blade.pro_t_C, Blade.pro_AoA, Blade.pro_cL, Blade.pro_cD] = read_pc_file(aerofoil_file);
        end
    end
end

