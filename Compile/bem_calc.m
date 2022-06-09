function [output, output_details, BEM] = bem_calc(file)
    
    [General, op_pts, Blade] = read_turbine_file(file);

    Blade.preflap = zeros(length(Blade.r),1);

    [output_details, output, BEM] = core_bem(General, op_pts, Blade);
    
    write_output(output_details, output, file);
    
end
