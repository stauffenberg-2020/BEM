function [output_details,output,BEM] = Turbine_caller(file)
%     file = 'G:\BEM\BEM\Data\NREL_5MW.txt';
    [General, op_pts, Blade] = read_turbine_file(file);

    Blade.preflap = zeros(length(Blade.r),1);

    AeroFlags.induction = General.induction; % 0 or 1, 0 = Induction Off, 1 = Induction On
    AeroFlags.tip_loss = General.tip_loss; % 0 or 1, 0 = Prandtl's tip loss correction Off, 1 = Prandtl's tip loss correction On
    AeroFlags.highCT = General.highCT; % 0 or 1 or 2, 0 = Off, 1 = As per HANSEN eqn. 6.38, 2 = As per HANSEN eqn. 6.37

    [output_details,output,BEM] = core_bem(General, op_pts, Blade);

    outfile = 'output.txt';
    fileID = fopen(outfile,'w');

    tmp_string = repmat('%s \t',1,size(output,2));
    tmp_unit = repmat('%s \t\t',1,size(output,2));
    tmp_float = repmat('%8.2f \t',1,size(output,2));
    tmp_string = strcat(tmp_string,'\n');
    tmp_unit = strcat(tmp_unit,'\n');
    tmp_float = strcat(tmp_float,'\n');
    
    fprintf(fileID,tmp_string,output_details{2,:});
    fprintf(fileID,tmp_unit,output_details{3,:});
    fprintf(fileID,tmp_float,output');
    fclose(fileID);
end
