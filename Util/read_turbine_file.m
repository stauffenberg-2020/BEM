function [General, op_pts, Blade, Control] = read_turbine_file(file)
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
    filebyline = cell(n_lines,1);
    
    for i=1:n_lines
        filebyline{i,1} = fgetl(fid);
    end
    fclose(fid);
    
    header_rows = zeros(n_lines,1);
    for i=1:n_lines
        header_rows(i,1) = strncmpi('#',filebyline{i},1);
    end
    header_row_ids = find(header_rows);
    no_h = 5; % Manually defined to have 5 headers, GENERAL, OPERATIONAL_SET-POINTS, BLADE_DETAILS, AEROFOIL_FILE, CONTROL
    header_titles = cell(1,no_h);
    for i=1:no_h
        header_titl = filebyline{header_row_ids(i)};
        header_titles{i} = extractAfter(header_titl,'# ');
        header_titles{i} = strtrim(header_titles{i});
    end
    
    % Pre-initilization of structures for MATLAB Coder compatibility
    General = struct('N',0,'rho',0,'induction',0,'tip_loss',0,'highCT',0);
    coder.cstructname(General,'General');
    Control = struct('CTR_flag',0,'Cp_file','','Ct_file','','Prat',0,'Vin',0,'Vrat',0,'Vout',0,'Oin',0,'Orat',0);
    coder.cstructname(Control,'Control');
    op_pts = struct('wsp',0,'pitch',0,'rpm',0);
    Blade = struct('r',0,'C',0,'t_C',0,'AeroTwist',0,'preflap',0,'pro_t_C',0,'pro_AoA',0,'pro_cL',0,'pro_cD',0);
    
    coder.varsize('op_pts.wsp','op_pts.pitch','op_pts.rpm');
    coder.varsize('Blade.r','Blade.C','Blade.t_C','Blade.AeroTwist','Blade.preflap','Blade.pro_t_C','Blade.pro_AoA','Blade.pro_cL','Blade.pro_cD');
        
    for i=1:length(header_titles)
        switch header_titles{i}
            case 'GENERAL'
                General.N = real(str2double(extractBefore((filebyline{header_row_ids(i)+1}),';')));
                General.rho = real(str2double(extractBefore((filebyline{header_row_ids(i)+2}),';')));
                General.induction = real(str2double(extractBefore((filebyline{header_row_ids(i)+3}),';')));
                General.tip_loss = real(str2double(extractBefore((filebyline{header_row_ids(i)+4}),';')));
                General.highCT = real(str2double(extractBefore((filebyline{header_row_ids(i)+5}),';')));
            case 'OPERATIONAL_SET-POINTS'
                headers = split_string(filebyline{header_row_ids(i)+1},' ');
                data = zeros(header_row_ids(i+1)-header_row_ids(i)-3,length(headers));
                for j=1:header_row_ids(i+1)-header_row_ids(i)-3
                    tline = split_string(filebyline{header_row_ids(i)+j-1+2},' ');
                    for k=1:length(headers)
                        data(j,k) = real(str2double(tline{k}));
                    end
                end
                op_pts.wsp=data(:,1);
                op_pts.pitch = data(:,2);
                op_pts.rpm = data(:,3);
            case 'BLADE_DETAILS'
                headers = split_string(filebyline{header_row_ids(i)+1},' ');
                data_bld = zeros(header_row_ids(i+1)-header_row_ids(i)-3,length(headers));
                for j=1:header_row_ids(i+1)-header_row_ids(i)-3
                    tline = split_string(filebyline{header_row_ids(i)+j-1+2},' ');
                    for k=1:length(headers)
                        data_bld(j,k) = real(str2double(tline{k}));
                    end
                end
                Blade.r=data_bld(:,1);
                Blade.C=data_bld(:,2);
                Blade.t_C=data_bld(:,3);
                Blade.AeroTwist=data_bld(:,4);
            case 'AEROFOIL_FILE'
                aerofoil_file = strtrim(filebyline{header_row_ids(i)+1});
                aerofoil_file = ['/home/runner/work/BEM/BEM',aerofoil_file];
                [Blade.pro_t_C, Blade.pro_AoA, Blade.pro_cL, Blade.pro_cD] = read_pc_file(aerofoil_file);
            case 'CONTROL'
                Control.CTR_flag = real(str2double(extractBefore((filebyline{header_row_ids(i)+1}),';')));
                Control.Cp_file = strtrim(extractBefore((filebyline{header_row_ids(i)+2}),';'));
                Control.Ct_file = strtrim(extractBefore((filebyline{header_row_ids(i)+3}),';'));
                Control.Prat = real(str2double(extractBefore((filebyline{header_row_ids(i)+4}),';')));
                Control.Vin = real(str2double(extractBefore((filebyline{header_row_ids(i)+5}),';')));
                Control.Vrat = real(str2double(extractBefore((filebyline{header_row_ids(i)+6}),';')));
                Control.Vout = real(str2double(extractBefore((filebyline{header_row_ids(i)+7}),';')));
                Control.Oin = real(str2double(extractBefore((filebyline{header_row_ids(i)+8}),';')));
                Control.Orat = real(str2double(extractBefore((filebyline{header_row_ids(i)+9}),';')));
        end
    end
end