function [output, output_details, BEM] = bem_calc(file)
    
    [General, op_pts, Blade, CTRL] = read_turbine_file(file);

    Blade.preflap = zeros(length(Blade.r),1);
    
    if CTRL.CTR_flag ==1
        [pitch_range, lambda_range, Cp] = read_Cp_Ct(CTRL.Cp_file);
        [op_pts.pitch, ~] = core_pitch_ctr(op_pts.wsp, op_pts.rpm, General.rho, Blade.r(end), pitch_range, lambda_range, Cp, CTRL.Oin, CTRL.Orat, CTRL.Prat, 0.944);
    end

    [output_details, output, BEM] = core_bem(General, op_pts, Blade);
    
    write_output(output_details, output, file);
    
end
