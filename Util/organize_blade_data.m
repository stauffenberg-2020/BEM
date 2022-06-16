function Blade = organize_blade_data(sec_r, sec_C, sec_t_C, sec_twist, pc)
    % This function organizes all the blade structural and aerodynamic data
    % in a structure for better data handling
    
    % r [m] is is the blade section table (m x 1 vector)
    % C [m] is the Chord corresponding to each of the sec_r (m x 1 vector)
    % t_C is the thickness to chord ratio corresponding to each of the sec_r (m x 1 vector)
    % AeroTwist [deg] is is the Aerodynamic twist at each of the sec_r (m x 1 vector)
    % pro_t_C is the K x 1 thickneses available
    % pro_AoA is the L x 1 angles of attack
    % pro_cL is L x K cL values corresponding to each angles of attack and t_C
    % pro_cD is L x K cD values corresponding to each angles of attack and t_C
    
    % Input data organising, Making sure values are in column vector
    r = sec_r(:);
    C = sec_C(:);
    t_C = sec_t_C(:);
    AeroTwist = sec_twist(:);
    
    % Organizing the blade data in a structure
    Blade.r = r;
    Blade.C = C;
    Blade.t_C = t_C;
    Blade.AeroTwist = AeroTwist;
    Blade.preflap = zeros(length(sec_r),1);
    
    [Blade.pro_t_C, Blade.pro_AoA, Blade.pro_cL, Blade.pro_cD] = read_pc_file(pc);
end